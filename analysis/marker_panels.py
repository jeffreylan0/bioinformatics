"""Curated retinal cell-population marker panels, shared by the figures.

Hand-picked canonical markers; no single-cell reference was available in this
repository, so these are literature markers rather than data-derived signatures
and each panel is small. They are used to read the direction of the composition
shift, not to deconvolve proportions.
"""
# Curated marker panels. Housekeeping is a negative control, not a cell type;
# reactive gliosis is a state rather than a population, and is included because
# it is the one unambiguously induced programme in this dataset.
PANELS = {
    "Reactive gliosis": ["Gfap", "Edn2", "Fgf2", "Lcn2", "Serpina3n", "C3",
                         "Stat3", "Socs3", "Atf3"],
    "Microglia": ["Cx3cr1", "Aif1", "P2ry12", "Csf1r", "Tmem119", "C1qa",
                  "C1qb", "Ctss"],
    "Rod": ["Rho", "Gnat1", "Pde6a", "Pde6b", "Pde6g", "Cnga1", "Cngb1", "Sag",
            "Rcvrn", "Grk1", "Rp1", "Prph2", "Rom1", "Nrl", "Nr2e3", "Rs1",
            "Prcd", "Rpgrip1", "Aipl1", "Gnb1", "Pdc", "Rgs9"],
    "Müller glia": ["Glul", "Rlbp1", "Slc1a3", "Apoe", "Clu", "Vim", "Aqp4",
                    "Dkk3", "Car14", "Sox9"],
    "Housekeeping": ["Actb", "Gapdh", "Tbp", "Hprt", "Ppia", "Rpl13a", "B2m",
                     "Ywhaz", "Sdha", "Pgk1"],
    "Amacrine": ["Gad1", "Gad2", "Slc6a9", "Tfap2a", "Tfap2b", "Chat",
                 "Slc32a1", "Pax6"],
    "Bipolar": ["Vsx2", "Grm6", "Trpm1", "Prkca", "Cabp5", "Grik1", "Isl1",
                "Vsx1", "Gnb3", "Nyx"],
    "Retinal ganglion": ["Rbpms", "Pou4f1", "Pou4f2", "Pou4f3", "Nefl", "Sncg",
                         "Thy1", "Nrn1", "Slc17a6"],
    "Horizontal": ["Lhx1", "Onecut1", "Onecut2", "Calb1"],
    "Cone": ["Opn1sw", "Opn1mw", "Arr3", "Pde6c", "Pde6h", "Gnat2", "Gngt2",
             "Cnga3", "Guca1a"],
}
