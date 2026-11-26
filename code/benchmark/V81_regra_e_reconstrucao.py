#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V81 - Fecha a regra de contagem e reconstroi os 26 com verificacao
      complexo a complexo.

Contexto: o V80 mostrou que a composicao da nuvem e "capsideo polimero,
sem o parceiro, sem HETATM" -- 11 de 15 residuos exatos na 6I2K, com os
outros 4 em exatamente -1 atomo.  Falta explicar esse -1.

FASE A  -  sobre os 290 residuos da 6I2K, nao 15.
  Grade de 8 regras de leitura/contagem:
     altloc {removida, mantida} x cos {cos(radians(45)), 0.7071} x raio {<=, <}
  mais a leitura direta do arquivo de montagem do RCSB (coordenadas
  arredondadas) com a melhor das 8.

CRITERIO DE ACEITACAO, fixado antes de ver o resultado:
  reproduz exatamente, OU diferenca maxima de 1 atomo E nenhum residuo
  trocando de lado no limiar de 106.
  Se nenhuma regra passar, o script PARA e nada e escrito.

FASE B  -  reconstroi os 26.  Cada complexo tem a propria pegada
  publicada no disco, entao cada um verifica a propria composicao em vez
  de herdar a da 6I2K.  Salva .npz com a mascara verificada registrada.

NAO sobrescreve nada.  Escreve so em /home/juan/V78_montagens/.
NAO classifica estrutura como invalida: o que nao fecha e reportado e
segue para decisao, nao descartado.

Uso:
    python3 V81_regra_e_reconstrucao.py
    python3 V81_regra_e_reconstrucao.py --so-fase-a
"""

import os
import sys
import glob
import math
import time
import numpy as np

sys.path.insert(0, "/home/juan")

try:
    import gemmi
except ImportError:
    sys.exit("gemmi nao instalado.   pip install gemmi")

try:
    from V78_reconstroi_montagens import monta, ocupacao_cone, CAMINHOS, SAIDA
except ImportError as e:
    sys.exit("nao achei o V78_reconstroi_montagens.py em /home/juan  (%s)" % e)

try:
    from scipy.spatial import cKDTree
except ImportError:
    sys.exit("scipy necessario.   pip install scipy")

SO_FASE_A = "--so-fase-a" in sys.argv

SEMIANG = 45.0
RAIO    = 30.0
LIMIAR  = 106            # limiar de ocupacao publicado
RCSB_6I2K = "/home/juan/V50_assembleias_rcsb/6i2k-assembly1.cif"

DIRS_PEGADA = [
    "/home/juan/repo_a71_V69_20261125T130138_06f6e831/data/per_residue",
    "/home/juan/V49_pegadas",
    "/home/juan/abertos",
    "/home/juan/entrega_virus_snca_ppdock/02_dados/pegadas",
]

AU_6I2K = "/home/juan/abertos/cif_virus/6i2k.cif"


def barra(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


# ------------------------------------------------------------- leitura ----

def le_estrutura(caminho, altloc=False):
    """altloc=False fica so com a conformacao A; True mantem todas."""
    st = gemmi.read_structure(caminho)
    st.setup_entities()
    if not altloc:
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
        "st": st,
        "xyz":  np.asarray(xyz, dtype=np.float64),
        "elem": np.asarray(elem),
        "cad":  np.asarray(cad),
        "sub":  np.asarray(sub),
        "seq":  np.asarray(seq, dtype=np.int32),
        "comp": np.asarray(comp),
        "nome": np.asarray(nome),
        "het":  np.asarray(het, dtype=bool),
    }


def acha_pegada(pid):
    for d in DIRS_PEGADA:
        if not os.path.isdir(d):
            continue
        for pad in ("%s_*_pegada_V49*.tsv", "%s_*_pegada*.tsv"):
            achados = sorted(glob.glob(os.path.join(d, pad % pid)))
            if achados:
                return achados[0]
    return None


def le_pegada(caminho):
    """Devolve (cadeia, {residuo: blocking})."""
    base = os.path.basename(caminho)
    partes = base.split("_")
    cadeia = partes[1] if len(partes) > 1 else None
    dados = {}
    with open(caminho, "r", errors="replace") as fh:
        cols = [c.strip().lower() for c in fh.readline().rstrip("\n").split("\t")]
        if "residue" not in cols or "blocking" not in cols:
            return cadeia, {}
        ir, ib = cols.index("residue"), cols.index("blocking")
        for linha in fh:
            p = linha.rstrip("\n").split("\t")
            if len(p) <= max(ir, ib):
                continue
            try:
                dados[int(float(p[ir]))] = int(float(p[ib]))
            except ValueError:
                continue
    return cadeia, dados


# ------------------------------------------------------------ contagem ----

def conta_cone(vert, nuvem, arv, centro, cosmin, inclusivo):
    """Mesma geometria do V78, com o gate angular e o teste de raio abertos
    como parametro, para testar as variantes de regra."""
    u = vert - centro
    u = u / np.linalg.norm(u)
    idx = arv.query_ball_point(vert, RAIO + 1e-9)
    if not idx:
        return 0
    d = nuvem[idx] - vert
    dist = np.sqrt((d * d).sum(axis=1))
    viz = dist > 1e-6
    cos = np.zeros(len(d))
    cos[viz] = (d[viz] @ u) / dist[viz]
    dentro = (dist <= RAIO) if inclusivo else (dist < RAIO)
    return int(np.count_nonzero(viz & dentro & (cos >= cosmin)))


def avalia(pos_ca, res_ca, ocup, nuvem, centro, cosmin, inclusivo):
    arv = cKDTree(nuvem)
    difs, flips = [], 0
    for p, r in zip(pos_ca, res_ca):
        r = int(r)
        if r not in ocup:
            continue
        esp = ocup[r]
        obt = conta_cone(p, nuvem, arv, centro, cosmin, inclusivo)
        difs.append(obt - esp)
        if (esp <= LIMIAR) != (obt <= LIMIAR):
            flips += 1
    d = np.asarray(difs, dtype=int)
    return {
        "n": len(d),
        "exatos": int((d == 0).sum()),
        "max": int(np.abs(d).max()) if len(d) else -1,
        "med": float(np.abs(d).mean()) if len(d) else float("nan"),
        "flips": flips,
        "neg": int((d < 0).sum()),
        "pos": int((d > 0).sum()),
    }


# ------------------------------------------------------- composicoes -----

def candidatos(mt, cad_ref):
    cad, het = mt["cad"], mt["het"]
    cadeias = sorted(set(cad.tolist()))
    polimero = [c for c in cadeias if het[cad == c].sum() < 0.5 * (cad == c).sum()]
    tam = sorted(polimero, key=lambda c: -(cad == c).sum())

    saida = [("sem HETATM", ~het),
             ("todas as cadeias", np.ones(len(cad), bool))]
    if len(tam) > 1 and tam[0] != cad_ref:
        saida.append(("sem HETATM, sem %s" % tam[0], ~het & (cad != tam[0])))
    if len(tam) > 2:
        fora = [c for c in tam[:2] if c != cad_ref]
        if fora:
            saida.append(("sem HETATM, sem %s" % ",".join(fora),
                          ~het & ~np.isin(cad, fora)))
    for c in tam:
        if c == cad_ref:
            continue
        saida.append(("sem HETATM, sem %s" % c, ~het & (cad != c)))
    if len(tam) >= 4:
        quatro = sorted(tam[:4])
        if cad_ref in quatro:
            saida.append(("sem HETATM, so %s" % ",".join(quatro),
                          ~het & np.isin(cad, quatro)))

    vistos, out = set(), []
    for nome, m in saida:
        if not m.any():
            continue
        k = (int(m.sum()), int(np.flatnonzero(m)[0]))
        if k in vistos:
            continue
        vistos.add(k)
        out.append((nome, m))
    return out, polimero


def aceita(r):
    """Criterio fixado antes de ver o resultado."""
    if r["n"] == 0:
        return False
    if r["exatos"] == r["n"]:
        return True
    return r["max"] <= 1 and r["flips"] == 0


# ---------------------------------------------------------- FASE A -------

def fase_a():
    barra("FASE A  -  a regra de contagem, sobre os 290 residuos da 6I2K")

    peg = acha_pegada("6I2K")
    if peg is None:
        print("  pegada da 6I2K nao encontrada.  PARANDO.")
        return None
    cad_ref, ocup = le_pegada(peg)
    print("  pegada : %s" % peg)
    print("  cadeia %s, %d residuos, blocking de %d a %d"
          % (cad_ref, len(ocup), min(ocup.values()), max(ocup.values())))
    perto = sum(1 for v in ocup.values() if abs(v - LIMIAR) <= 2)
    print("  residuos a menos de 2 atomos do limiar %d: %d" % (LIMIAR, perto))

    resultados = []
    for altloc in (False, True):
        au = le_estrutura(AU_6I2K, altloc=altloc)
        mt, err = monta(au)
        if mt is None:
            print("  reconstrucao falhou (altloc=%s): %s" % (altloc, err))
            continue
        masc = (~mt["het"]) & (mt["cad"] != "E")
        nuvem = mt["xyz"][masc]
        centro = nuvem.mean(axis=0)
        ref = (mt["op"] == max(mt["id_op"], 0)) & (mt["cad"] == cad_ref) & (mt["nome"] == "CA")
        pos_ca, res_ca = mt["xyz"][ref], mt["seq"][ref]

        # o V78 e a referencia: conferir que a versao parametrizada coincide
        if not altloc:
            arv = cKDTree(nuvem)
            for p in pos_ca[:5]:
                a = ocupacao_cone(p, nuvem, centro)
                b = conta_cone(p, nuvem, arv, centro, math.cos(math.radians(SEMIANG)), True)
                if a != b:
                    print("  AVISO: versao parametrizada diverge do V78 (%d vs %d)" % (a, b))
                    break

        for nome_cos, cosmin in (("cos(45)", math.cos(math.radians(SEMIANG))),
                                 ("0.7071 ", 0.7071)):
            for nome_r, inclusivo in (("<=", True), ("< ", False)):
                r = avalia(pos_ca, res_ca, ocup, nuvem, centro, cosmin, inclusivo)
                r["rotulo"] = "altloc %s | %s | dist %s" % (
                    "mantida " if altloc else "removida", nome_cos, nome_r)
                r["altloc"] = altloc
                r["cosmin"] = cosmin
                r["inclusivo"] = inclusivo
                resultados.append(r)

    print("\n  %-36s %7s %6s %6s %7s %6s %6s"
          % ("regra", "exatos", "de", "max", "media", "neg", "flips"))
    for r in resultados:
        print("  %-36s %7d %6d %6d %7.3f %6d %6d"
              % (r["rotulo"], r["exatos"], r["n"], r["max"], r["med"],
                 r["neg"], r["flips"]))

    # leitura direta do arquivo do RCSB, com a melhor das 8
    melhor = max(resultados, key=lambda r: (r["exatos"], -r["med"]))
    if os.path.exists(RCSB_6I2K):
        print("\n  mesma regra, lendo o arquivo de montagem do RCSB:")
        rc = le_estrutura(RCSB_6I2K, altloc=melhor["altloc"])
        # no arquivo do RCSB as copias ja estao expandidas; a cadeia de
        # referencia e a que mantem o nome original
        masc = (~rc["het"]) & (~np.char.startswith(rc["cad"].astype(str), "E"))
        nuvem = rc["xyz"][masc]
        centro = nuvem.mean(axis=0)
        ref = (rc["cad"] == cad_ref) & (rc["nome"] == "CA")
        if ref.sum() == 0:
            print("    cadeia %s nao encontrada com esse nome no arquivo do RCSB" % cad_ref)
        else:
            r = avalia(rc["xyz"][ref], rc["seq"][ref], ocup, nuvem, centro,
                       melhor["cosmin"], melhor["inclusivo"])
            print("    exatos %d de %d   max %d   media %.3f   flips %d"
                  % (r["exatos"], r["n"], r["max"], r["med"], r["flips"]))

    vencedora = None
    for r in sorted(resultados, key=lambda x: (-x["exatos"], x["med"])):
        if aceita(r):
            vencedora = r
            break

    if vencedora is None:
        barra("NENHUMA regra passou no criterio.  Nada foi escrito.")
        print("  criterio: exato, ou max 1 atomo E zero trocas no limiar %d." % LIMIAR)
        print("  melhor foi '%s': %d de %d exatos, max %d, %d trocas."
              % (melhor["rotulo"], melhor["exatos"], melhor["n"],
                 melhor["max"], melhor["flips"]))
        return None

    print("\n  PASSOU: %s" % vencedora["rotulo"])
    print("  %d de %d exatos, max %d atomo, %d trocas no limiar %d."
          % (vencedora["exatos"], vencedora["n"], vencedora["max"],
             vencedora["flips"], LIMIAR))
    return vencedora


# ---------------------------------------------------------- FASE B -------

def fase_b(regra):
    barra("FASE B  -  reconstroi os 26, cada um verificado pela propria pegada")
    os.makedirs(SAIDA, exist_ok=True)
    print("  regra: %s" % regra["rotulo"])
    print("  saida: %s\n" % SAIDA)

    todos = dict(CAMINHOS)
    todos["6I2K"] = AU_6I2K

    print("  %-6s %9s %-30s %9s %6s %6s %7s"
          % ("pdb", "atomos", "composicao", "exatos/de", "max", "flips", "s"))
    relato = []
    for pid, caminho in sorted(todos.items()):
        t0 = time.time()
        if not os.path.exists(caminho):
            print("  %-6s ARQUIVO AUSENTE: %s" % (pid, caminho))
            relato.append((pid, "arquivo ausente", None))
            continue
        peg = acha_pegada(pid)
        if peg is None:
            print("  %-6s sem pegada publicada, nao da para verificar" % pid)
            relato.append((pid, "sem pegada", None))
            continue
        cad_ref, ocup = le_pegada(peg)
        if not ocup:
            print("  %-6s pegada sem coluna blocking" % pid)
            relato.append((pid, "pegada ilegivel", None))
            continue

        try:
            au = le_estrutura(caminho, altloc=regra["altloc"])
            mt, err = monta(au)
        except Exception as e:                                  # noqa: BLE001
            print("  %-6s falhou na leitura: %s" % (pid, e))
            relato.append((pid, "erro: %s" % e, None))
            continue
        if mt is None:
            print("  %-6s %s" % (pid, err))
            relato.append((pid, err, None))
            continue

        ref_op = max(mt["id_op"], 0)
        if cad_ref not in set(mt["cad"].tolist()):
            print("  %-6s cadeia %s da pegada nao existe na estrutura" % (pid, cad_ref))
            relato.append((pid, "cadeia %s ausente" % cad_ref, None))
            continue
        ref = (mt["op"] == ref_op) & (mt["cad"] == cad_ref) & (mt["nome"] == "CA")
        pos_ca, res_ca = mt["xyz"][ref], mt["seq"][ref]

        cands, _ = candidatos(mt, cad_ref)
        melhor_r, melhor_nome, melhor_masc = None, None, None
        for nome_c, masc in cands:
            nuvem = mt["xyz"][masc]
            if len(nuvem) < 1000:
                continue
            r = avalia(pos_ca, res_ca, ocup, nuvem, nuvem.mean(axis=0),
                       regra["cosmin"], regra["inclusivo"])
            if melhor_r is None or (r["exatos"], -r["med"]) > (melhor_r["exatos"], -melhor_r["med"]):
                melhor_r, melhor_nome, melhor_masc = r, nome_c, masc
            if r["exatos"] == r["n"]:
                break

        dt = time.time() - t0
        ok = aceita(melhor_r) if melhor_r else False
        print("  %-6s %9d %-30s %4d/%-4d %6d %6d %7.0f%s"
              % (pid, len(mt["xyz"]), melhor_nome[:30], melhor_r["exatos"],
                 melhor_r["n"], melhor_r["max"], melhor_r["flips"], dt,
                 "" if ok else "   <-- NAO PASSOU"))
        relato.append((pid, melhor_nome, melhor_r))

        if not ok:
            continue

        nuvem_masc = melhor_masc
        destino = os.path.join(SAIDA, "%s_montagem.npz" % pid)
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
            mascara=nuvem_masc,
            cad_ref=np.bytes_(str(cad_ref).encode()),
            id_op=np.int32(mt["id_op"]),
            n_op=np.int32(mt["n_op"]),
            centro=mt["xyz"][nuvem_masc].mean(axis=0),
            composicao=np.bytes_(melhor_nome.encode()),
            regra=np.bytes_(regra["rotulo"].encode()),
            fonte=np.bytes_(caminho.encode()),
            pegada=np.bytes_(peg.encode()),
        )

    barra("RESUMO")
    passou = [r for r in relato if r[2] is not None and aceita(r[2])]
    print("  verificados e salvos : %d de %d" % (len(passou), len(todos)))
    falhou = [r for r in relato if not (r[2] is not None and aceita(r[2]))]
    if falhou:
        print("  nao passaram:")
        for pid, nome, r in falhou:
            if r is None:
                print("    %-6s %s" % (pid, nome))
            else:
                print("    %-6s %s  (%d/%d exatos, max %d, %d trocas)"
                      % (pid, nome, r["exatos"], r["n"], r["max"], r["flips"]))
        print("\n  Isto nao os torna invalidos.  E a minha reconstrucao que nao")
        print("  reproduziu a contagem publicada neles; a decisao e sua.")


def main():
    barra("V81 - regra de contagem e reconstrucao verificada")
    print("  semi-angulo %.0f graus   raio %.0f A   limiar %d" % (SEMIANG, RAIO, LIMIAR))
    regra = fase_a()
    if regra is None:
        return
    if SO_FASE_A:
        barra("--so-fase-a: fase B nao executada.")
        return
    fase_b(regra)


if __name__ == "__main__":
    main()
