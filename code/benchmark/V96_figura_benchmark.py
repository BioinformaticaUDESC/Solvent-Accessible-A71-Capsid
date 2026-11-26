#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V96 - Figura do benchmark: forest plot das diferencas pareadas.

Le /home/juan/V83_benchmark/V94_tabela4.tsv, que e a fonte canonica
gravada pelo V94.  Nenhum numero digitado aqui.

Um grafico, um arquivo.  Sem titulo dentro da figura, sem blocos de
texto, rotulos de eixo em ingles.  PNG a 600 dpi e TIFF a 300 dpi com
compressao LZW, largura de 180 mm.

Nao sobrescreve nada: escreve em /home/juan/V96_figura_benchmark/.

Uso:
    python3 V96_figura_benchmark.py
"""

import os
import sys
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FONTE = "/home/juan/V83_benchmark/V94_tabela4.tsv"
SAIDA = "/home/juan/V96_figura_benchmark"
LARGURA_MM = 180.0

ROTULOS = {
    "Accessible area, 1.4 A probe":  "Accessible area, 1.4 Å probe",
    "Accessible area, 5.0 A probe":  "Accessible area, 5.0 Å probe",
    "Accessible area, 10.0 A probe": "Accessible area, 10.0 Å probe",
    "Residue depth":                 "Residue depth",
    "Radial percentile":             "Radial percentile",
    "Half-sphere exposure":          "Half-sphere exposure",
}


def main():
    if not os.path.exists(FONTE):
        sys.exit("nao achei %s.  rode o V94 antes." % FONTE)
    os.makedirs(SAIDA, exist_ok=True)

    linhas = []
    with open(FONTE) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not r.get("dif_pareada"):
                continue
            linhas.append((r["descritor"],
                           float(r["dif_pareada"]),
                           float(r["dif_lo"]),
                           float(r["dif_hi"]),
                           r["cone_vence"], r["de"]))
    if not linhas:
        sys.exit("o TSV nao tem linhas com diferenca pareada.")

    # maior diferenca em cima
    linhas.sort(key=lambda x: x[1], reverse=True)
    print("  %d descritores lidos de %s" % (len(linhas), FONTE))
    for d, m, lo, hi, v, de in linhas:
        print("    %-32s %+0.4f [%+0.4f, %+0.4f]  %s/%s" % (d, m, lo, hi, v, de))

    pol = LARGURA_MM / 25.4
    alt = 0.52 * len(linhas) + 1.25
    fig, ax = plt.subplots(figsize=(pol, alt))

    y = list(range(len(linhas)))[::-1]
    for yi, (d, m, lo, hi, v, de) in zip(y, linhas):
        ax.plot([lo, hi], [yi, yi], color="0.25", linewidth=1.4,
                solid_capstyle="butt", zorder=2)
        ax.plot([lo, lo], [yi - 0.12, yi + 0.12], color="0.25", linewidth=1.4, zorder=2)
        ax.plot([hi, hi], [yi - 0.12, yi + 0.12], color="0.25", linewidth=1.4, zorder=2)
        ax.plot([m], [yi], marker="o", markersize=5.0, color="0.10",
                zorder=3, linestyle="none")
        ax.text(hi + 0.006, yi, "%s/%s" % (v, de), va="center", ha="left",
                fontsize=8, color="0.25")

    ax.axvline(0.0, color="0.55", linewidth=0.9, linestyle="--", zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([ROTULOS.get(d, d) for d, *_ in linhas], fontsize=9)
    ax.set_xlabel("Paired difference in AUC (outward-cone occupancy − descriptor)",
                  fontsize=9)
    ax.tick_params(axis="x", labelsize=8)
    ax.set_ylim(-0.6, len(linhas) - 0.4)
    xmax = max(h for _, _, _, h, _, _ in linhas)
    ax.set_xlim(min(0.0, min(l for _, _, l, _, _, _ in linhas)) - 0.012,
                xmax + 0.045)
    for lado in ("top", "right", "left"):
        ax.spines[lado].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.tight_layout()

    png = os.path.join(SAIDA, "Figure_5_benchmark_paired_differences.png")
    tif = os.path.join(SAIDA, "Figure_5_benchmark_paired_differences.tif")
    fig.savefig(png, dpi=600, facecolor="white")
    fig.savefig(tif, dpi=300, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

    print("\n  %s  (%.1f MB)" % (png, os.path.getsize(png) / 1e6))
    print("  %s  (%.1f MB)" % (tif, os.path.getsize(tif) / 1e6))
    print("\n  Legenda sugerida, para o texto e nao para dentro da figura:")
    print("    Figure 5. Paired difference in AUC between outward-cone occupancy")
    print("    (45 deg, 30 A) and each alternative descriptor, computed within")
    print("    each of 22 assembled complexes. Points are the mean paired")
    print("    difference and bars the 95% interval from 5000 resamples of")
    print("    complexes. The dashed line marks no difference. Numbers at the")
    print("    right give the complexes in which occupancy ranked higher.")


if __name__ == "__main__":
    main()
