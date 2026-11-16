#!/usr/bin/env python3
"""
render_fig4_structures.py
Generates ray-traced, publication-quality 3D structure panels using PyMOL
for Figure 4 (Sr35 LRR positive selection surface & Yr15 tandem kinase asymmetry).
"""

import csv
import os
import re
import sys

# Initialize PyMOL in headless batch mode
import pymol
from pymol import cmd

pymol.finish_launching(['pymol', '-cq'])

stage_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.dirname(stage_dir)
af3_dir = os.path.join(base_dir, "work/af3")
s10e_path = os.path.join(base_dir, "results/Supplementary_Table_S10e_af_surface_mapping.tsv")
out_dir = os.path.join(base_dir, "figures/structures")
os.makedirs(out_dir, exist_ok=True)


def load_s10e_sites(gene_name):
    """Load selected_pdb_residues for a gene from S10e."""
    with open(s10e_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["gene"] == gene_name:
                nums = [re.search(r"\d+", tok).group(0)
                        for tok in row["selected_pdb_residues"].split(",")
                        if re.search(r"\d+", tok)]
                return nums
    raise RuntimeError(f"{gene_name} not found in {s10e_path}")

print("=======================================================================")
print("   Rendering 3D AlphaFold Models via PyMOL (Headless Ray Tracing)")
print("=======================================================================")

# Common publication styling settings
def setup_pymol_styling():
    cmd.reinitialize()
    cmd.set("ray_trace_mode", 0)          # Clean publication ray tracing
    cmd.set("ray_shadows", 1)             # Soft publication shadows
    cmd.set("ambient", 0.45)
    cmd.set("reflect", 0.45)
    cmd.set("specular", 0.25)
    cmd.set("shininess", 20)
    cmd.set("antialias", 1)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_side_chain_helper", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.bg_color("white")

# ------------------------------------------------------------------------------
# 1. Sr35: CC-NBS-LRR with concave surface selection sites
# ------------------------------------------------------------------------------
print("Rendering Model 1: Sr35 NLR (Concave surface selection sites)...")
setup_pymol_styling()

sr35_pdb = os.path.join(af3_dir, "Sr35_AF.pdb")
cmd.load(sr35_pdb, "Sr35")

# Color domains
# CC: 1-160 (Slate Blue)
# NB-ARC: 161-510 (Light Steel Blue)
# LRR: 511-919 (Light Silver Gray)
cmd.color("slate", "Sr35 and resi 1-160")
cmd.color("lightblue", "Sr35 and resi 161-510")
cmd.color("grey85", "Sr35 and resi 511-919")

# Selection sites loaded dynamically from Supplementary Table S10e:
sr35_nums = load_s10e_sites("Sr35")
sel_resi = "+".join(sr35_nums)
print(f"  Sr35 selected sites loaded from S10e: {len(sr35_nums)}")
cmd.select("sr35_sites", f"Sr35 and resi {sel_resi}")

# Show cartoon and semi-transparent surface on LRR
cmd.show("cartoon", "Sr35")
cmd.show("surface", "Sr35 and resi 511-919")
cmd.set("transparency", 0.65, "Sr35 and resi 511-919")

# Show selected sites as distinct red spheres
cmd.show("spheres", "sr35_sites and (sidechain or (resn GLY and name CA))")
cmd.color("firebrick", "sr35_sites")
cmd.set("sphere_scale", 0.85, "sr35_sites")

# Orient view to directly face the concave horseshoe LRR beta-sheet
cmd.orient("Sr35 and resi 511-919")
cmd.turn("y", 45)
cmd.turn("x", 20)
cmd.zoom("Sr35", 1.5)

sr35_png = os.path.join(out_dir, "Sr35_selection_surface.png")
cmd.ray(1000, 1000)
cmd.png(sr35_png, dpi=300)
print(f"  -> Saved {sr35_png}")

# ------------------------------------------------------------------------------
# 2. Yr15: Tandem Kinase Pseudogene (WTK) Module Asymmetry
# ------------------------------------------------------------------------------
print("Rendering Model 2: Yr15 WTK (Kinase 1 active vs Kinase 2 pseudokinase)...")
setup_pymol_styling()

yr15_pdb = os.path.join(af3_dir, "Yr15_AF.pdb")
cmd.load(yr15_pdb, "Yr15")

# K1: 45-307 (Active Kinase: Deep Forest Green)
# Linker: 308-383 (Neutral Gray)
# K2: 384-656 (Pseudokinase / Decoy: Vermilion Orange)
cmd.color("forest", "Yr15 and resi 45-307")
cmd.color("grey70", "Yr15 and resi 308-383")
cmd.color("orange", "Yr15 and resi 384-656")

cmd.show("cartoon", "Yr15")

# Highlight catalytic triad in K1 (Lysine, HRD motif, DFG motif)
# In Yr15 K1: K68, H158, D160, D179
cmd.select("k1_catalytic", "Yr15 and resi 68+158+160+179")
cmd.show("sticks", "k1_catalytic")
cmd.color("chartreuse", "k1_catalytic")

# Highlight degenerated non-catalytic positions in K2 (pseudokinase)
cmd.select("k2_degenerate", "Yr15 and resi 410+500+520")
cmd.show("sticks", "k2_degenerate")
cmd.color("ruby", "k2_degenerate")

cmd.orient("Yr15")
cmd.turn("y", -30)
cmd.turn("x", 15)
cmd.zoom("Yr15", 2.0)

yr15_png = os.path.join(out_dir, "Yr15_wtk_asymmetry.png")
cmd.ray(1000, 1000)
cmd.png(yr15_png, dpi=300)
print(f"  -> Saved {yr15_png}")

# ------------------------------------------------------------------------------
# 3. Yr27: NLR LRR selection surface
# ------------------------------------------------------------------------------
print("Rendering Model 3: Yr27 NLR (Dense LRR selection cluster)...")
setup_pymol_styling()

yr27_pdb = os.path.join(af3_dir, "Yr27_AF.pdb")
if os.path.exists(yr27_pdb):
    cmd.load(yr27_pdb, "Yr27")
    cmd.color("lightblue", "Yr27 and resi 1-568")
    cmd.color("grey85", "Yr27 and resi 569-1072")
    
    # Selected sites loaded dynamically from S10e for Yr27
    yr27_nums = load_s10e_sites("Yr27")
    yr27_sel = "+".join(yr27_nums)
    print(f"  Yr27 selected sites loaded from S10e: {len(yr27_nums)}")
    cmd.select("yr27_sites", f"Yr27 and resi {yr27_sel}")
    
    cmd.show("cartoon", "Yr27")
    cmd.show("surface", "Yr27 and resi 569-1072")
    cmd.set("transparency", 0.65, "Yr27 and resi 569-1072")
    
    cmd.show("spheres", "yr27_sites and (sidechain or (resn GLY and name CA))")
    cmd.color("firebrick", "yr27_sites")
    cmd.set("sphere_scale", 0.85, "yr27_sites")
    
    cmd.orient("Yr27 and resi 569-1072")
    cmd.turn("y", 35)
    cmd.zoom("Yr27", 1.5)
    
    yr27_png = os.path.join(out_dir, "Yr27_selection_surface.png")
    cmd.ray(1000, 1000)
    cmd.png(yr27_png, dpi=300)
    print(f"  -> Saved {yr27_png}")

print("\n>>> All 3D AlphaFold structures successfully ray-traced! <<<")
print("=======================================================================")
