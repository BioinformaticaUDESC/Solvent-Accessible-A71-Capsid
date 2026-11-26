#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V78 - Reconstroi as 25 montagens que faltam e VERIFICA a reconstrucao
      contra os numeros de ocupacao no cone ja publicados.

NAO sobrescreve nada.  Escreve so em /home/juan/V78_montagens/.
NAO executa simulacao: expande os operadores depositados no proprio mmCIF.
NAO classifica estrutura como invalida; se algo nao bate, ele reporta e para.

Ordem:
  1. 6I2K primeiro.  Reconstroi a partir da unidade assimetrica e compara
     (a) com a montagem do RCSB que ja esta salva em V50_assembleias_rcsb/
     (b) com a ocupacao no cone da tabela por residuo publicada.
     Testa variantes de composicao (todas as cadeias, e dropando uma de
     cada vez) para descobrir QUAL regra reproduz o numero publicado,
     em vez de assumir uma.
  2. So se (b) bater, reconstroi os 25 que faltam e salva .npz.

Uso:
    python3 V78_reconstroi_montagens.py            # verifica e, se bater, reconstroi tudo
    python3 V78_reconstroi_montagens.py --so-teste # para depois da verificacao
"""

import os
import sys
import math
import glob
import time
import numpy as np

try:
    import gemmi
except ImportError:
    sys.exit("gemmi nao instalado.   pip install gemmi")

SAIDA       = "/home/juan/V78_montagens"
RCSB_6I2K   = "/home/juan/V50_assembleias_rcsb/6i2k-assembly1.cif"
SEMIANG     = 45.0      # semi-angulo do cone, graus
RAIO        = 30.0      # alcance da busca, A
SO_TESTE    = "--so-teste" in sys.argv

# caminhos exatos reportados pelo inventario do V77 -- nao sao chute
CAMINHOS = {
    "1V9U": "/home/juan/abertos/cif_virus/1v9u.cif",
    "3DPR": "/home/juan/abertos/cif_virus/3dpr.cif",
    "3J8F": "/home/juan/abertos/cif_virus/3j8f.cif",
    "6ADL": "/home/juan/cif_proveniencia/6ADL.cif",
    "6ADM": "/home/juan/cif_proveniencia/6ADM.cif",
    "6ADR": "/home/juan/abertos/cif_virus/6adr.cif",
    "6CX1": "/home/juan/abertos/cif_virus/6cx1.cif",
    "6EIT": "/home/juan/abertos/cif_virus/6eit.cif",
    "6ILJ": "/home/juan/cif_proveniencia/6ILJ.cif",
    "6ILK": "/home/juan/cif_proveniencia/6ILK.cif",
    "6PPO": "/home/juan/cif_proveniencia/6PPO.cif",
    "6PSF": "/home/juan/abertos/cif_virus/6psf.cif",
    "6SNW": "/home/juan/abertos/cif_virus/6snw.cif",
    "7BZT": "/home/juan/abertos/cif_virus/7bzt.cif",
    "7BZU": "/home/juan/abertos/cif_virus/7bzu.cif",
    "7C9W": "/home/juan/cif_proveniencia/7C9W.cif",
    "7VXZ": "/home/juan/abertos/cif_virus/7vxz.cif",
    "7VY5": "/home/juan/abertos/cif_virus/7vy5.cif",
    "7VY6": "/home/juan/cif_proveniencia/7VY6.cif",
    "7VYK": "/home/juan/abertos/cif_virus/7vyk.cif",
    "7VYL": "/home/juan/abertos/cif_virus/7vyl.cif",
    "7VYM": "/home/juan/abertos/cif_virus/7vym.cif",
    "7W14": "/home/juan/cif_proveniencia/7W14.cif",
    "7XXA": "/home/juan/abertos/cif_virus/7xxa.cif",
    "8B8R": "/home/juan/cif_proveniencia/8B8R.cif",
}


def barra(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


# ---------------------------------------------------------------- leitura ---

def acha_au_6i2k():
    for d in ("/home/juan/abertos/cif_virus", "/home/juan/cif_proveniencia"):
        for n in ("6i2k.cif", "6I2K.cif"):
            p = os.path.join(d, n)
            if os.path.exists(p):
                return p
    achados = glob.glob("/home/juan/**/6i2k.cif", recursive=True)
    achados += glob.glob("/home/juan/**/6I2K.cif", recursive=True)
    return achados[0] if achados else None


def le_au(caminho):
    """Le a unidade assimetrica.  Remove agua, hidrogenio e confs alternativas."""
    st = gemmi.read_structure(caminho)
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.remove_waters()
    model = st[0]

    xyz, elem, cad, sub, seq, comp, nome, het = [], [], [], [], [], [], [], []
    for chain in model:
        for res in chain:
            h = (res.het_flag == "H")
            for at in res:
                xyz.append((at.pos.x, at.pos.y, at.pos.z))
                elem.append(at.element.name)
                cad.append(chain.name)
                sub.append(res.subchain)
                seq.append(res.seqid.num)
                comp.append(res.name)
                nome.append(at.name)
                het.append(h)

    return {
        "st":   st,
        "xyz":  np.asarray(xyz, dtype=np.float64),
        "elem": np.asarray(elem),
        "cad":  np.asarray(cad),
        "sub":  np.asarray(sub),
        "seq":  np.asarray(seq, dtype=np.int32),
        "comp": np.asarray(comp),
        "nome": np.asarray(nome),
        "het":  np.asarray(het, dtype=bool),
    }


def _transf(t):
    try:
        M = np.asarray(t.mat.tolist(), dtype=np.float64)
    except AttributeError:
        M = np.asarray([[t.mat[i][j] for j in range(3)] for i in range(3)], dtype=np.float64)
    try:
        v = np.asarray(t.vec.tolist(), dtype=np.float64)
    except AttributeError:
        v = np.asarray([t.vec.x, t.vec.y, t.vec.z], dtype=np.float64)
    return M, v


def monta(au):
    """Expande a unidade assimetrica pelos operadores depositados.

    Devolve arrays da particula inteira com, para cada atomo, o indice do
    operador que o gerou.  copia 0 = primeiro operador do primeiro gerador.
    """
    st = au["st"]
    if not st.assemblies:
        return None, "nenhuma assembly depositada no arquivo"

    asm = st.assemblies[0]

    blocos_xyz, blocos_op = [], []
    blocos_idx = []
    iop = 0
    for gen in asm.generators:
        chs = list(gen.chains)
        sbs = list(gen.subchains)
        if sbs:
            masc = np.isin(au["sub"], sbs)
        elif chs:
            masc = np.isin(au["cad"], chs)
        else:
            masc = np.ones(len(au["xyz"]), dtype=bool)

        if not masc.any():
            continue

        base = au["xyz"][masc]
        idx  = np.nonzero(masc)[0]
        for op in gen.operators:
            M, v = _transf(op.transform)
            blocos_xyz.append(base @ M.T + v)
            blocos_op.append(np.full(len(base), iop, dtype=np.int32))
            blocos_idx.append(idx)
            iop += 1

    if not blocos_xyz:
        return None, "os geradores da assembly nao selecionaram nenhum atomo"

    xyz = np.concatenate(blocos_xyz, axis=0)
    opi = np.concatenate(blocos_op,  axis=0)
    idx = np.concatenate(blocos_idx, axis=0)

    # identidade: qual operador reproduz a unidade assimetrica como depositada
    id_op = -1
    k = 0
    for gen in asm.generators:
        for op in gen.operators:
            M, v = _transf(op.transform)
            if np.allclose(M, np.eye(3), atol=1e-6) and np.allclose(v, 0.0, atol=1e-6):
                id_op = k
                break
            k += 1
        if id_op >= 0:
            break

    return {
        "xyz":   xyz,
        "op":    opi,
        "orig":  idx,                 # indice do atomo correspondente na AU
        "elem":  au["elem"][idx],
        "cad":   au["cad"][idx],
        "sub":   au["sub"][idx],
        "seq":   au["seq"][idx],
        "comp":  au["comp"][idx],
        "nome":  au["nome"][idx],
        "het":   au["het"][idx],
        "n_op":  iop,
        "id_op": id_op,
    }, None


# -------------------------------------------------------------- geometria ---

def ocupacao_cone(vert, nuvem, centro, semiang=SEMIANG, raio=RAIO):
    """Atomos da nuvem dentro do cone de vertice `vert`, eixo radial para fora.

    Exclui apenas o atomo do proprio vertice (teste de distancia minima);
    cadeia lateral propria e vizinhos sequenciais CONTAM, como no publicado.
    """
    u = vert - centro
    n = np.linalg.norm(u)
    if n < 1e-9:
        return -1
    u = u / n
    d = nuvem - vert
    dist = np.sqrt((d * d).sum(axis=1))
    viz = dist > 1e-6
    cos = np.zeros(len(d))
    cos[viz] = (d[viz] @ u) / dist[viz]
    lim = math.cos(math.radians(semiang))
    return int(np.count_nonzero(viz & (dist <= raio) & (cos >= lim)))


# ----------------------------------------------- tabela por residuo 6I2K ---

ALVOS_OCUP = ("ocup", "occup", "cone", "cono", "atoms_cone", "n_cone", "bloq")
ALVOS_RES  = ("residue", "resid", "resseq", "res_num", "residuo", "seqid")


def acha_tabela_6i2k():
    """Procura a tabela por residuo da 6I2K e devolve {residuo: ocupacao}."""
    cands = []
    for ext in ("tsv", "csv", "txt"):
        cands += glob.glob("/home/juan/**/*6i2k*.%s" % ext, recursive=True)
        cands += glob.glob("/home/juan/**/*6I2K*.%s" % ext, recursive=True)
    cands = sorted(set(cands))

    melhor = None
    for c in cands:
        try:
            with open(c, "r", errors="replace") as fh:
                cab = fh.readline().rstrip("\n")
        except OSError:
            continue
        sep = "\t" if "\t" in cab else ","
        cols = [x.strip().lower() for x in cab.split(sep)]
        col_r = next((i for i, c2 in enumerate(cols) if c2 in ALVOS_RES), None)
        col_o = next((i for i, c2 in enumerate(cols)
                      if any(a in c2 for a in ALVOS_OCUP)), None)
        if col_r is None or col_o is None:
            continue
        dados = {}
        with open(c, "r", errors="replace") as fh:
            fh.readline()
            for linha in fh:
                p = linha.rstrip("\n").split(sep)
                if len(p) <= max(col_r, col_o):
                    continue
                try:
                    dados[int(float(p[col_r]))] = int(float(p[col_o]))
                except ValueError:
                    continue
        if dados and (melhor is None or len(dados) > len(melhor[2])):
            melhor = (c, cols[col_o], dados)
    return melhor


# ------------------------------------------------------------ verificacao ---

def verifica_6i2k():
    barra("1.  VERIFICACAO  -  6I2K, a unica com montagem ja salva")

    au_path = acha_au_6i2k()
    if au_path is None:
        print("  unidade assimetrica da 6I2K nao encontrada.  nao da para verificar.")
        return False
    print("  unidade assimetrica : %s" % au_path)

    au = le_au(au_path)
    print("  atomos na AU        : %d   cadeias: %s"
          % (len(au["xyz"]), ", ".join(sorted(set(au["cad"])))))

    print("\n  inventario de cadeias da AU")
    print("    %-6s %8s %9s %10s %12s" % ("cadeia", "atomos", "residuos", "het", "faixa resid"))
    for c in sorted(set(au["cad"])):
        m = au["cad"] == c
        rr = au["seq"][m]
        print("    %-6s %8d %9d %10d %6d-%-6d"
              % (c, m.sum(), len(set(rr.tolist())), au["het"][m].sum(), rr.min(), rr.max()))

    mt, err = monta(au)
    if mt is None:
        print("\n  reconstrucao falhou: %s" % err)
        return False
    print("\n  operadores depositados : %d" % mt["n_op"])
    print("  operador identidade    : %s"
          % ("nenhum (a AU nao aparece intacta na montagem)" if mt["id_op"] < 0
             else "indice %d" % mt["id_op"]))
    print("  atomos na reconstrucao : %d" % len(mt["xyz"]))

    # (a) contra a montagem do RCSB
    if os.path.exists(RCSB_6I2K):
        print("\n  (a) contra a montagem do RCSB")
        t0 = time.time()
        rc = le_au(RCSB_6I2K)
        print("      atomos no RCSB     : %d" % len(rc["xyz"]))
        print("      atomos reconstruido: %d" % len(mt["xyz"]))
        c1 = mt["xyz"].mean(axis=0)
        c2 = rc["xyz"].mean(axis=0)
        print("      centroides distam  : %.3f A" % float(np.linalg.norm(c1 - c2)))
        try:
            from scipy.spatial import cKDTree
            arv = cKDTree(rc["xyz"])
            dd, _ = arv.query(mt["xyz"], k=1)
            print("      vizinho mais proximo: mediana %.4f A   p99 %.4f A   max %.4f A"
                  % (float(np.median(dd)), float(np.percentile(dd, 99)), float(dd.max())))
        except ImportError:
            print("      scipy ausente, pulando a comparacao atomo a atomo")
        print("      (%.1f s)" % (time.time() - t0))
    else:
        print("\n  (a) montagem do RCSB nao encontrada em %s, pulando" % RCSB_6I2K)

    # (b) contra a ocupacao publicada
    print("\n  (b) contra a ocupacao no cone da tabela por residuo")
    tab = acha_tabela_6i2k()
    if tab is None:
        print("      nenhuma tabela por residuo da 6I2K encontrada sob /home/juan.")
        print("      sem ela nao da para confirmar que a reconstrucao e a mesma")
        print("      que gerou os numeros publicados.  PARANDO.")
        return False
    arq, colname, ocup = tab
    print("      tabela : %s" % arq)
    print("      coluna : %s   (%d residuos)" % (colname, len(ocup)))

    # cadeia de referencia: a da copia identidade cujos residuos mais
    # coincidem com os da tabela
    ref_op = mt["id_op"] if mt["id_op"] >= 0 else 0
    m_ref = mt["op"] == ref_op
    melhor_cad, melhor_n = None, -1
    for c in sorted(set(mt["cad"][m_ref])):
        mm = m_ref & (mt["cad"] == c)
        comum = len(set(mt["seq"][mm].tolist()) & set(ocup.keys()))
        if comum > melhor_n:
            melhor_cad, melhor_n = c, comum
    print("      cadeia de referencia : %s   (%d residuos em comum com a tabela)"
          % (melhor_cad, melhor_n))
    if melhor_n < 10:
        print("      coincidencia baixa demais para verificar.  PARANDO.")
        return False

    # vertices: CA da cadeia de referencia que estao na tabela
    m_ca = m_ref & (mt["cad"] == melhor_cad) & (mt["nome"] == "CA")
    res_ca = mt["seq"][m_ca]
    pos_ca = mt["xyz"][m_ca]
    sel = [i for i, r in enumerate(res_ca) if int(r) in ocup]
    if len(sel) > 15:
        passo = len(sel) // 15
        sel = sel[::passo][:15]
    print("      testando %d residuos" % len(sel))

    # variantes de composicao
    cadeias = sorted(set(mt["cad"]))
    variantes = [("todas as cadeias", np.ones(len(mt["xyz"]), dtype=bool))]
    if not mt["het"].all() and mt["het"].any():
        variantes.append(("sem HETATM", ~mt["het"]))
    for c in cadeias:
        if len(cadeias) > 1:
            variantes.append(("sem a cadeia %s" % c, mt["cad"] != c))

    print()
    print("      %-24s %8s %8s %10s" % ("variante", "iguais", "testados", "erro med"))
    vencedora = None
    for nome_v, masc in variantes:
        nuvem = mt["xyz"][masc]
        centro = nuvem.mean(axis=0)
        iguais, erros = 0, []
        for i in sel:
            esperado = ocup[int(res_ca[i])]
            obtido = ocupacao_cone(pos_ca[i], nuvem, centro)
            erros.append(abs(obtido - esperado))
            if obtido == esperado:
                iguais += 1
        med = float(np.mean(erros)) if erros else float("nan")
        print("      %-24s %8d %8d %10.1f" % (nome_v, iguais, len(sel), med))
        if iguais == len(sel) and vencedora is None:
            vencedora = (nome_v, masc)

    if vencedora is None:
        print("\n      NENHUMA variante reproduziu a ocupacao publicada.")
        print("      A diferenca e de composicao da nuvem (quais atomos entram na")
        print("      contagem), nao de geometria.  Isso precisa sair do")
        print("      filtro_externo.py, nao de palpite.  PARANDO.")
        print("      Me manda:  grep -n -i 'centro\\|centroid\\|cone\\|0.7071\\|45' \\")
        print("                   /home/juan/abertos/filtro_externo.py | head -40")
        return False

    print("\n      reproduz o publicado com: %s" % vencedora[0])
    print("      a reconstrucao e equivalente a que gerou os numeros do paper.")
    return vencedora[0]


# -------------------------------------------------------------- producao ---

def reconstroi_todos(regra):
    barra("2.  RECONSTRUCAO  -  os 25 que faltam")
    os.makedirs(SAIDA, exist_ok=True)
    print("  saida: %s" % SAIDA)
    print("  regra de composicao confirmada na 6I2K: %s" % regra)
    print()

    tempos = []
    for i, (pid, caminho) in enumerate(sorted(CAMINHOS.items()), 1):
        destino = os.path.join(SAIDA, "%s_montagem.npz" % pid)
        if os.path.exists(destino):
            print("  %2d/25  %-6s ja existe, pulando" % (i, pid))
            continue
        if not os.path.exists(caminho):
            print("  %2d/25  %-6s ARQUIVO AUSENTE: %s" % (i, pid, caminho))
            continue

        t0 = time.time()
        try:
            au = le_au(caminho)
            mt, err = monta(au)
        except Exception as e:                                  # noqa: BLE001
            print("  %2d/25  %-6s FALHOU na leitura: %s" % (i, pid, e))
            continue
        if mt is None:
            print("  %2d/25  %-6s FALHOU: %s" % (i, pid, err))
            continue

        np.savez_compressed(
            destino,
            xyz=mt["xyz"].astype(np.float32),
            op=mt["op"],
            elem=np.char.encode(mt["elem"], "ascii"),
            cad=np.char.encode(mt["cad"], "ascii"),
            seq=mt["seq"],
            comp=np.char.encode(mt["comp"], "ascii"),
            nome=np.char.encode(mt["nome"], "ascii"),
            het=mt["het"],
            id_op=np.int32(mt["id_op"]),
            n_op=np.int32(mt["n_op"]),
            centro=mt["xyz"].mean(axis=0),
            fonte=np.bytes_(caminho.encode()),
        )
        dt = time.time() - t0
        tempos.append(dt)
        mb = os.path.getsize(destino) / 1e6
        print("  %2d/25  %-6s  %7d atomos   %2d operadores   %5.1f MB   %5.1f s"
              % (i, pid, len(mt["xyz"]), mt["n_op"], mb, dt))

        if len(tempos) == 1:
            resta = (len(CAMINHOS) - 1) * dt
            print("         >>> ao ritmo do primeiro, os 25 levam ~%.0f min" % (resta / 60.0))

    if tempos:
        print("\n  total %.1f min   mediana %.1f s por complexo"
              % (sum(tempos) / 60.0, float(np.median(tempos))))


def main():
    barra("V78 - reconstrucao verificada das montagens")
    print("  semi-angulo %.0f graus   raio %.0f A   (os mesmos do publicado)" % (SEMIANG, RAIO))

    regra = verifica_6i2k()
    if not regra:
        barra("PAROU NA VERIFICACAO.  Nada foi escrito.")
        return
    if SO_TESTE:
        barra("--so-teste: verificacao passou, reconstrucao nao foi executada.")
        return
    reconstroi_todos(regra)
    barra("pronto.  o benchmark (V79) le direto de %s" % SAIDA)


if __name__ == "__main__":
    main()
