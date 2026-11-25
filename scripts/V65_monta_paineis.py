#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V65 — monta as figuras individuais nos quatro compostos de submissao.

    cd /home/juan && python3 V65_monta_paineis.py

Entra em /home/juan/V58_figuras, sai em /home/juan/V65_submissao.
Nao altera nenhuma figura individual: so le e compoe.

Proporcao de cada imagem preservada (nenhuma distorcao).  Painel alinhado
pelo topo, altura comum por linha, espacamento uniforme, fundo branco puro.
Letra de painel em negrito no canto superior esquerdo.

Largura: 180 mm, pagina inteira da Wiley.
Saida: PNG 600 dpi e TIFF 300 dpi com compressao LZW.

Os nomes sao procurados com ou sem sufixo de download (_2, _3), e o script
imprime qual arquivo usou em cada painel.  Se faltar um, ele diz e monta a
figura com os que existem.
"""

import os, sys, glob, re

try:
    from PIL import Image, ImageDraw, ImageFont
except ImportError:
    sys.exit("Pillow ausente: pip install --user Pillow")

BASE    = "/home/juan"
ENTRADA = f"{BASE}/V58_figuras"
SAIDA   = f"{BASE}/V65_submissao"
os.makedirs(SAIDA, exist_ok=True)

DPI_PNG, DPI_TIFF = 600, 300
LARG_MM = 180.0
W = int(round(LARG_MM / 25.4 * DPI_PNG))      # 4252 px
MARGEM, GAP = 36, 64
ALT_LETRA = 86                                 # faixa reservada a letra

FIGURAS = [
    ("Figure_1_method_and_curation",
     [("Fig0_cone_definition", "A"),
      ("fluxograma_curadoria",  "B")]),
    ("Figure_2_calibration_and_evaluation",
     [("Fig1_occupancy_distribution_calibration", "A"),
      ("Fig2_threshold_sweep_heldout",            "B"),
      ("Fig3_occupancy_vs_radial_percentile",     "C")]),
    ("Figure_3_sasa_paradox_and_site22",
     [("Fig4_sasa_retention_vs_radial",           "A"),
      ("Fig5_site22_occupancy_across_deposits",   "B")]),
    ("Figure_4_structural_context_3VBS",
     [("Fig8_protomer_site22",        "A"),
      ("Fig6_capsid_outer_surface",   "B"),
      ("Fig7_capsid_equatorial_section", "C")]),
]

# nomes alternativos aceitos para o fluxograma, que nao sei como voce salvou
ALIASES = {"fluxograma_curadoria": ["fluxograma*", "*curadoria*", "*curation*",
                                    "*flowchart*", "*Figure_02*", "*fluxo*"]}


def acha(base):
    """arquivo mais recente que casa a base, com ou sem sufixo _2 _3 ..."""
    pats = [f"{base}.png", f"{base}_*.png", f"{base}.jpg", f"{base}_*.jpg",
            f"{base}.tif", f"{base}_*.tif"]
    for a in ALIASES.get(base, []):
        pats += [f"{a}.png", f"{a}.jpg"]
    achados = []
    for p in pats:
        achados += glob.glob(os.path.join(ENTRADA, p))
        achados += glob.glob(os.path.join(ENTRADA, "**", p), recursive=True)
    achados = [a for a in set(achados) if os.path.isfile(a)]
    if not achados:
        return None
    # prefere nome exato; senao o mais recente
    exato = [a for a in achados
             if re.fullmatch(rf"{re.escape(base)}\.(png|jpg|tif)",
                             os.path.basename(a), re.I)]
    return exato[0] if exato else max(achados, key=os.path.getmtime)


def carrega(p):
    im = Image.open(p)
    if im.mode in ("RGBA", "LA", "P"):
        fundo = Image.new("RGB", im.size, "white")
        im = im.convert("RGBA")
        fundo.paste(im, mask=im.split()[-1])
        im = fundo
    return im.convert("RGB")


def fonte(px):
    import matplotlib.font_manager as fm
    for nome in ("DejaVu Sans:bold", "Arial:bold", "Helvetica:bold"):
        try:
            cam = fm.findfont(fm.FontProperties(
                family=nome.split(":")[0], weight="bold"), fallback_to_default=False)
            return ImageFont.truetype(cam, px)
        except Exception:
            continue
    try:
        return ImageFont.truetype(fm.findfont(fm.FontProperties(weight="bold")), px)
    except Exception:
        return ImageFont.load_default()


FT = fonte(ALT_LETRA - 20)

print("=" * 74)
print("V65 — montagem dos compostos")
print("=" * 74)
print(f"\nentrada: {ENTRADA}")
print(f"saida:   {SAIDA}")
print(f"largura: {LARG_MM:.0f} mm = {W} px a {DPI_PNG} dpi\n")

for nome_saida, paineis in FIGURAS:
    print("-" * 74)
    print(nome_saida)
    itens = []
    for base, letra in paineis:
        p = acha(base)
        if p is None:
            print(f"  [{letra}] AUSENTE  {base}")
            continue
        im = carrega(p)
        itens.append((im, letra, p))
        print(f"  [{letra}] {im.size[0]}x{im.size[1]}  {os.path.basename(p)}")

    if not itens:
        print("  nenhum painel disponivel; figura nao montada")
        continue

    # painel unico perde a letra, como voce pediu
    unico = len(itens) == 1
    if unico:
        print("  painel unico: montada sem letra de painel")

    n = len(itens)
    disp = W - 2 * MARGEM - (n - 1) * GAP
    soma = sum(im.size[0] / im.size[1] for im, _, _ in itens)
    H = int(round(disp / soma))                       # altura comum da linha

    larguras = []
    for im, _, _ in itens:
        larguras.append(int(round(H * im.size[0] / im.size[1])))
    # ajuste de arredondamento na ultima
    larguras[-1] += disp - sum(larguras)

    topo = 0 if unico else ALT_LETRA
    lona = Image.new("RGB", (W, H + topo + 2 * MARGEM), "white")
    dz = ImageDraw.Draw(lona)

    x = MARGEM
    for (im, letra, _), w in zip(itens, larguras):
        lona.paste(im.resize((w, H), Image.LANCZOS), (x, MARGEM + topo))
        if not unico:
            dz.text((x, MARGEM), letra, font=FT, fill="black")
        x += w + GAP

    png = os.path.join(SAIDA, f"{nome_saida}.png")
    lona.save(png, "PNG", dpi=(DPI_PNG, DPI_PNG))
    print(f"  -> {os.path.basename(png)}  {lona.size[0]}x{lona.size[1]}  "
          f"{os.path.getsize(png)/1e6:.1f} MB")

    fator = DPI_TIFF / DPI_PNG
    pequena = lona.resize((int(lona.size[0]*fator), int(lona.size[1]*fator)),
                          Image.LANCZOS)
    tif = os.path.join(SAIDA, f"{nome_saida}.tif")
    pequena.save(tif, "TIFF", compression="tiff_lzw",
                 dpi=(DPI_TIFF, DPI_TIFF))
    print(f"  -> {os.path.basename(tif)}  {pequena.size[0]}x{pequena.size[1]}  "
          f"{os.path.getsize(tif)/1e6:.1f} MB")

print("\n" + "=" * 74)
for p in sorted(glob.glob(os.path.join(SAIDA, "*"))):
    print(f"  {os.path.getsize(p)/1e6:7.1f} MB  {os.path.basename(p)}")
print(f"\ncompostos em {SAIDA}")
print("as figuras individuais nao foram alteradas")
print("=" * 74)
