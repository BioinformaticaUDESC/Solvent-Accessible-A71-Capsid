#!/usr/bin/env python3
"""Perfil radial de VP1 no capsídeo montado: quais resíduos ficam na face externa.

Para cada resíduo de VP1 calcula, no capsídeo inteiro:
  - raio do CA em relação ao centro do capsídeo
  - posição na espessura da casca (0% interna, 100% externa)
  - direção da cadeia lateral (para fora / para dentro)
  - átomos que o cobrem num cone radial externo

Reporta os sítios candidatos da reconstrução ancestral e, como controle, a
alça GH de VP1, que é onde SCARB2 liga.

    pip install --user gemmi numpy matplotlib
    python3 perfil_radial_vp1.py 3VBS.pdb
    python3 perfil_radial_vp1.py 3VBS.pdb --sitios 22,98,225,282,289
"""
import argparse
import sys
from pathlib import Path

import numpy as np

try:
    import gemmi
except ImportError:
    sys.exit("falta gemmi:  pip install --user gemmi")

AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
      "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
BACKBONE = {"N", "CA", "C", "O", "OXT"}
GH = range(208, 223)          # alça GH de VP1: onde SCARB2 liga
CONE_COS = np.cos(np.radians(45.0))
RAIO = 30.0
EXTERNO = 70.0                # limiar de "face externa", em % da casca


def do_cabecalho(pdb):
    mol, atual = {}, None
    for linha in pdb.read_text(errors="ignore").splitlines():
        if not linha.startswith("COMPND"):
            if linha.startswith(("ATOM", "HETATM", "SEQRES")):
                break
            continue
        corpo = linha[10:].strip().rstrip(";")
        if corpo.upper().startswith("MOLECULE:"):
            atual = corpo.split(":", 1)[1].strip()
        elif corpo.upper().startswith("CHAIN:") and atual:
            for c in corpo.split(":", 1)[1].split(","):
                mol[c.strip()] = atual
    return mol


def acha_vp1(st, pdb):
    cab = do_cabecalho(pdb)
    tam = {ch.name: sum(1 for r in ch if r.name in AA) for ch in st[0]}
    for c in sorted(tam):
        print(f"  {c}: {tam[c]:4d} res  {cab.get(c, '(sem COMPND)')}")
    alvos = [c for c in sorted(cab) if "VP1" in cab[c].upper()]
    if alvos:
        return alvos[0]
    cand = [c for c, n in tam.items() if abs(n - 297) <= 10]
    if len(cand) == 1:
        return cand[0]
    sys.exit("nao consegui identificar VP1; verifique o PDB")


def expande(st):
    try:
        if st.assemblies:
            m = gemmi.make_assembly(st.assemblies[0], st[0],
                                    gemmi.HowToNameCopiedChain.AddNumber)
            if sum(1 for _ in m) > sum(1 for _ in st[0]):
                return m, "BIOMT"
    except Exception:
        pass
    ops = [o for o in st.ncs if not o.given]
    if not ops:
        return None, None
    g = gemmi.Model("1")
    for ch in st[0]:
        g.add_chain(ch.clone())
    for i, op in enumerate(ops, 1):
        for ch in st[0]:
            novo = ch.clone()
            novo.name = f"{ch.name}~{i}"
            for r in novo:
                for a in r:
                    a.pos = gemmi.Position(op.tr.apply(a.pos))
            g.add_chain(novo)
    return g, f"MTRIX ({len(ops)} ops)"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb")
    ap.add_argument("--sitios", default="22,98,225,282,289",
                    help="sitios candidatos da ASR, separados por virgula")
    ap.add_argument("--out", default=None, help="prefixo dos arquivos de saida")
    ap.add_argument("--limiar", type=float, default=EXTERNO,
                    help="%% da casca para contar como externo (padrao 70)")
    args = ap.parse_args()

    pdb = Path(args.pdb)
    if not pdb.is_file():
        sys.exit(f"nao encontrei {pdb}")
    sitios = [int(x) for x in args.sitios.split(",") if x.strip()]
    pref = args.out or pdb.stem
    limiar = args.limiar

    st = gemmi.read_structure(str(pdb))
    st.setup_entities()
    print("cadeias:")
    vp1 = acha_vp1(st, pdb)
    print(f"VP1 = cadeia {vp1}\n")

    grande, origem = expande(st)
    if grande is None:
        sys.exit("sem operadores de simetria; nao da para montar o capsideo")
    n_ch = sum(1 for _ in grande)
    print(f"capsideo montado por {origem}: {n_ch} cadeias")

    coords = np.array([[a.pos.x, a.pos.y, a.pos.z]
                       for ch in grande for r in ch for a in r], dtype=np.float32)
    centro = coords.mean(axis=0)
    print(f"centro: {centro.round(1)}  sobre {len(coords)} atomos")

    ca_r = np.array([np.linalg.norm([a.pos.x, a.pos.y, a.pos.z] - centro)
                     for ch in grande for r in ch if r.name in AA
                     for a in r if a.name == "CA"], dtype=np.float32)
    r_int, r_ext = float(ca_r.min()), float(ca_r.max())
    print(f"casca: {r_int:.1f} a {r_ext:.1f} A (espessura {r_ext - r_int:.1f} A)\n")
    ca_sort = np.sort(ca_r)

    # copia identidade de VP1 no capsideo
    prim = next(r for r in st[0][vp1] if r.name in AA)
    p0 = np.array([prim[0].pos.x, prim[0].pos.y, prim[0].pos.z])
    alvo_ch, d_min = None, 1e9
    for ch in grande:
        for r in ch:
            if r.name in AA:
                d = np.linalg.norm([r[0].pos.x, r[0].pos.y, r[0].pos.z] - p0)
                if d < d_min:
                    alvo_ch, d_min = ch, d
                break
    print(f"copia identidade: cadeia {alvo_ch.name} (desvio {d_min:.3f} A)\n")

    linhas = []
    print("calculando perfil radial...", end="", flush=True)
    for res in alvo_ch:
        if res.name not in AA:
            continue
        ca = next((a for a in res if a.name == "CA"), res[0])
        pos = np.array([ca.pos.x, ca.pos.y, ca.pos.z])
        r_ca = float(np.linalg.norm(pos - centro))
        casca = 100.0 * float(np.searchsorted(ca_sort, r_ca)) / len(ca_sort)

        lado = [a for a in res if a.name not in BACKBONE]
        if lado:
            c = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in lado]).mean(axis=0)
            delta = float(np.linalg.norm(c - centro)) - r_ca
        else:
            delta = 0.0

        u = (pos - centro) / r_ca
        d = coords - pos
        dist = np.linalg.norm(d, axis=1)
        m = (dist > 1e-6) & (dist < RAIO)
        proj = (d[m] @ u) / dist[m]
        bloq = int((proj > CONE_COS).sum())

        bfac = float(np.mean([a.b_iso for a in res])) if len(res) else 0.0
        linhas.append((res.seqid.num, res.name, r_ca, casca, delta, bloq, bfac))
    print(f" {len(linhas)} residuos")

    bs = np.sort(np.array([l[6] for l in linhas]))
    def b_pct(v):
        return 100.0 * float(np.searchsorted(bs, v)) / len(bs)

    # Linha de base para ordem: a alca GH, sitio de ligacao do receptor
    # validado experimentalmente nesta mesma estrutura. Comparar B de residuo
    # exposto contra a cadeia inteira e errado: o nucleo enterrado e muito
    # mais ordenado, e o criterio rejeitaria a propria alca GH.
    b_gh = float(np.mean([l[6] for l in linhas if l[0] in GH])) or 1.0
    TOL = 1.5

    def viavel(c, b, bf):
        """Exposto e nao mais movel que o sitio de ligacao conhecido."""
        return c >= limiar and b < 50 and bf <= TOL * b_gh

    tsv = Path(f"{pref}_perfil_radial_VP1.tsv")
    with open(tsv, "w") as fh:
        fh.write("residue\tresname\tradius_A\tshell_pct\tsidechain_delta_A\t"
                 "blocking_atoms\tmean_bfactor\tbfactor_pct\t"
                 "b_over_GHloop\texposed\tviable\n")
        for n, nome, r, c, dl, b, bf in linhas:
            fh.write(f"{n}\t{nome}\t{r:.2f}\t{c:.2f}\t{dl:+.2f}\t{b}\t"
                     f"{bf:.1f}\t{b_pct(bf):.1f}\t{bf / b_gh:.2f}\t"
                     f"{'yes' if c >= limiar and b < 50 else 'no'}\t"
                     f"{'yes' if viavel(c, b, bf) else 'no'}\n")
    print(f"gravado {tsv}")
    print(f"B-factor da cadeia: mediana {np.median(bs):.1f}, "
          f"p75 {np.percentile(bs, 75):.1f}, maximo {bs[-1]:.1f}")
    print(f"B-factor da alca GH (referencia de ordem): {b_gh:.1f}; "
          f"limite de viabilidade {TOL:.1f}x = {TOL * b_gh:.1f}")

    por_num = {n: (nome, r, c, dl, b, bf) for n, nome, r, c, dl, b, bf in linhas}
    gh = [por_num[n][2] for n in GH if n in por_num]
    gh_b = [por_num[n][5] for n in GH if n in por_num]

    print(f"\n{'sitio':<7}{'res':<5}{'casca':>7}{'lado':>8}{'cob':>6}"
          f"{'B':>8}{'B/GH':>7}  veredito")
    print("-" * 74)
    for n in sitios:
        if n not in por_num:
            print(f"{n:<7}{'-':<5}{'nao modelado':>30}")
            continue
        nome, r, c, dl, b, bf = por_num[n]
        direc = "fora" if dl > 0 else "dentro"
        if viavel(c, b, bf):
            ver = "VIAVEL"
        elif c >= limiar and b < 50:
            ver = "exposto, mal ordenado"
        else:
            ver = "nao acessivel"
        print(f"{n:<7}{nome:<5}{c:6.1f}%{direc:>8}{b:6d}{bf:8.1f}"
              f"{bf / b_gh:7.2f}  {ver}")
    if gh:
        print("-" * 74)
        print(f"{'GH':<7}{'208-222':<5}{np.mean(gh):5.1f}%{'':>8}{'':>6}"
              f"{np.mean(gh_b):8.1f}{1.00:7.2f}  controle positivo (referencia)")

    # -------------------------------- o filtro PREVE o sitio de ligacao?
    # Atencao ao desenho: o eixo de ordem foi calibrado na propria alca GH,
    # entao usa-lo para "prever" a alca GH seria circular. A validacao abaixo
    # usa SOMENTE exposicao, que nao depende da alca, e trata os residuos da
    # alca GH como positivos conhecidos.
    pos = [l for l in linhas if l[0] in GH]
    neg = [l for l in linhas if l[0] not in GH]
    if pos and neg:
        maior = sum(1 for p in pos for q in neg if p[3] > q[3])
        iguais = sum(1 for p in pos for q in neg if p[3] == q[3])
        auc = (maior + 0.5 * iguais) / (len(pos) * len(neg))
        ordenado = sorted(linhas, key=lambda l: -l[3])
        posicoes = [i + 1 for i, l in enumerate(ordenado) if l[0] in GH]
        n10 = max(1, len(linhas) // 10)
        no_top = sum(1 for p in posicoes if p <= n10)
        esperado = len(pos) * n10 / len(linhas)
        print(f"\nO filtro preve a alca GH? (apenas o eixo de exposicao)")
        print(f"  AUC                       {auc:.3f}   (0,5 = acaso)")
        print(f"  posicoes da alca GH       {min(posicoes)}-{max(posicoes)} "
              f"de {len(linhas)}, mediana {int(np.median(posicoes))}")
        print(f"  no decil superior         {no_top} de {len(pos)} "
              f"(esperado por acaso {esperado:.1f})")
        if auc >= 0.8:
            print("  -> exposicao sozinha ja discrimina o sitio de ligacao")
        elif auc >= 0.65:
            print("  -> discriminacao moderada")
        else:
            print("  -> exposicao sozinha NAO discrimina; o filtro serve para")
            print("     excluir sitios impossiveis, nao para prever sitios de ligacao")
        print("  (o eixo de ordem nao entra aqui: foi calibrado nesta mesma alca)")

    # ------------------------------------------------------------ figuras
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib ausente; so o TSV foi gravado)")
        return

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                         "savefig.dpi": 300, "pdf.fonttype": 42})
    nums = [l[0] for l in linhas]
    cascas = [l[3] for l in linhas]
    bfs = [l[6] for l in linhas]

    fig, ax = plt.subplots(figsize=(6.4, 2.6))
    ax.axhspan(limiar, 100, color="#2B6CB0", alpha=0.07, zorder=0)
    ax.axvspan(GH.start, GH.stop - 1, color="#2F855A", alpha=0.12, zorder=0)
    ax.plot(nums, cascas, lw=0.8, color="#444444", zorder=2)
    for n in sitios:
        if n in por_num:
            ax.plot(n, por_num[n][2], "o", ms=5, color="#C05621",
                    mec="white", mew=0.8, zorder=4)
            ax.annotate(str(n), (n, por_num[n][2]), textcoords="offset points",
                        xytext=(0, 7), ha="center", fontsize=6.8,
                        color="#C05621")
    ax.set_xlabel("VP1 residue", fontsize=8)
    ax.set_ylabel("Position in shell (%)", fontsize=8)
    ax.set_xlim(min(nums) - 5, max(nums) + 5)
    ax.set_ylim(-3, 112)
    ax.set_yticks([0, 50, 100])
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ax.spines.values():
        s.set_linewidth(0.6)
    ax.tick_params(length=2.5, width=0.6, labelsize=7.5)
    ax.yaxis.grid(True, color="#DCDCDC", lw=0.5)
    ax.set_axisbelow(True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{pref}_Figure_radial_profile.{ext}",
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 2.2))
    ax.axvspan(GH.start, GH.stop - 1, color="#2F855A", alpha=0.12, zorder=0)
    ax.plot(nums, bfs, lw=0.8, color="#444444", zorder=2)
    ax.axhline(b_gh, color="#2F855A", lw=0.7, ls="--", zorder=1)
    ax.axhline(TOL * b_gh, color="#767676", lw=0.6, ls=":", zorder=1)
    for n in sitios:
        if n in por_num:
            ax.plot(n, por_num[n][5], "o", ms=5, color="#C05621",
                    mec="white", mew=0.8, zorder=4)
            ax.annotate(str(n), (n, por_num[n][5]), textcoords="offset points",
                        xytext=(0, 7), ha="center", fontsize=6.8,
                        color="#C05621")
    ax.set_xlabel("VP1 residue", fontsize=8)
    ax.set_ylabel("Mean B-factor", fontsize=8)
    ax.set_xlim(min(nums) - 5, max(nums) + 5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)
    ax.tick_params(length=2.5, width=0.6, labelsize=7.5)
    ax.yaxis.grid(True, color="#DCDCDC", lw=0.5)
    ax.set_axisbelow(True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{pref}_Figure_bfactor_profile.{ext}",
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"\ngravado {pref}_Figure_radial_profile.png / .pdf")
    print(f"gravado {pref}_Figure_bfactor_profile.png / .pdf")


if __name__ == "__main__":
    main()
