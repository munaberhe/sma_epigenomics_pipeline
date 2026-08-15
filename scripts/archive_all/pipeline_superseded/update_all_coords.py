with open('/data/home/bt25018/sma_epigenomics_pipeline/scripts/pipeline/20_master_locus_plots.R', 'r') as f:
    content = f.read()

coords = {
    "SEMA3C":   ("chr7",  80799100,  80829100,  80765501,  80815949),
    "SIGMAR1":  ("chr9",  34625116,  34655116,  34634739,  34637844),
    "RELL2":    ("chr5",  141600122, 141630122, 141589123, 141689422),
    "DDIT4L":   ("chr4",  100157963, 100187963, 100141264, 100241563),
    "MRPS2":    ("chr9",  135477616, 135507616, 135453567, 135553866),
    "GNG14":    ("chr19", 12676918,  12706918,  12641619,  12742218),
    "RNA5S13":  ("chr1",  228662668, 228692668, 228587769, 228688068),
    "KIAA1656": ("chr22", 30349684,  30379684,  30331935,  30432234),
    "TCEAL4":   ("chrX",  103579302, 103609302, 103536053, 103636352),
    "IRF8":     ("chr16", 85879034,  85909034,  85865935,  85966234),
    "USP27X":   ("chrX",  49909302,  49939302,  49832753,  49933052),
    "USP7":     ("chr16", 8912084,   8942084,   8840635,   8940934),
    "KDM1A":    ("chr1",  23034568,  23064568,  23035569,  23135868),
    "PAX5_ASO": ("chr9",  37046116,  37076116,  36992067,  37092366),
    "PAX5_VPA": ("chr9",  37046116,  37076116,  36992067,  37092366),
    "CAMK2A":   ("chr5",  150245072, 150275072, 150164623, 150264922),
    "EPHB1":    ("chr3",  135030969, 135060969, 135041420, 135141719),
    "ZDHHC22":  ("chr14", 77081022,  77111022,  77092523,  77192822),
}

for gene, (chr_, zs, ze, gs, ge) in coords.items():
    # find current start/end and replace with zoom coords
    import re
    pattern = f'(list\\(name="{gene}", chr="{chr_}", start=)(\\d+)(, end=)(\\d+)(,)'
    replacement = f'\\g<1>{zs}\\g<3>{ze}\\g<5>'
    new_content = re.sub(pattern, replacement, content)
    if new_content != content:
        content = new_content
        print(f"Updated {gene}: {zs}-{ze}")
    else:
        print(f"NOT FOUND: {gene}")

# Update output directories to new structure
content = content.replace(
    'OUT_ANNOTATED <- "results/dmr/plots/annotated"',
    'OUT_CANDIDATES <- "results/thesis_figures/locus_candidates"\nOUT_SMN2_ALL   <- "results/thesis_figures/locus_smn2"'
)
content = content.replace(
    'dir.create(OUT_ANNOTATED,', 
    'dir.create(OUT_CANDIDATES,\ndir.create(OUT_SMN2_ALL,'
)

with open('/data/home/bt25018/sma_epigenomics_pipeline/scripts/pipeline/20_master_locus_plots.R', 'w') as f:
    f.write(content)
print("Done")
