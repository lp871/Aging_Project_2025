## ============================================================
## METHODS (for manuscript)
## ============================================================
## Gene set enrichment analysis of aging-associated transcriptional changes.
## To interrogate the biological processes underlying retinal aging, we assembled
## a curated, species-matched gene-set catalog organized into nine aging-related
## functional categories: (C1) genome and epigenome maintenance, (C2) proteostasis,
## autophagy and the lysosome, (C3) mitochondrial metabolism and nutrient sensing,
## (C4) cellular stress, senescence and cell death, (C5) inflammation and innate
## immunity, (C6) intercellular communication and the glial niche, (C7) retinal cell
## identity and specialized function, (C8) regeneration and developmental plasticity,
## and (C9) extracellular matrix, vascular function and tissue architecture. For each
## category we manually selected representative MSigDB Hallmark gene sets and Gene
## Ontology Biological Process (GO:BP) terms. Gene sets were retrieved for the
## corresponding species with the msigdbr package (Hallmark collection "H" and C5
## GO:BP subcollection); GO terms were mapped to GO identifiers using GO.db and
## AnnotationDbi. We further incorporated two curated cellular senescence signatures
## (SenMayo and the Reactome cellular senescence set). Gene sets were compiled into a
## unified catalog providing both term-level pathways and category-level gene unions,
## which was used as input for gene set enrichment analysis (fgsea) on genes ranked
## by their age-associated expression change.
## ============================================================

#######
####### we will caluculate Aging GSEA #######
#######
####### first we need to create our own GO terms dataset #########
#######

1. Genome/epigenome
2. Proteostasis/autophagy
3. Metabolism/mitochondria
4. Stress/senescence
5. Inflammation
6. Intercellular communication
7. Retinal function
8. Regeneration/plasticity
9. ECM/vascular remodeling


#######
####### The nine functional categories provide a unified biological framework ########
#######

conda activate clusterProfiler
R

library(msigdbr)
library(babelgene)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tibble)
library(readr)
library(GO.db)
library(AnnotationDbi)
library(fgsea)

########

## ============================================================
## 1. Define nine aging functional categories
## ============================================================

aging_category_order <- c(
  "C1_Genome_epigenome",
  "C2_Proteostasis_autophagy",
  "C3_Metabolism_mitochondria",
  "C4_Stress_senescence_cell_death",
  "C5_Inflammation_immunity",
  "C6_Intercellular_communication",
  "C7_Retinal_cell_function",
  "C8_Regeneration_plasticity",
  "C9_ECM_vascular_architecture"
)

aging_category_labels <- tribble(
  ~category, ~category_label,

  "C1_Genome_epigenome",
  "Genome and epigenome maintenance",

  "C2_Proteostasis_autophagy",
  "Proteostasis, autophagy and lysosome",

  "C3_Metabolism_mitochondria",
  "Mitochondrial metabolism and nutrient sensing",

  "C4_Stress_senescence_cell_death",
  "Cellular stress, senescence and cell death",

  "C5_Inflammation_immunity",
  "Inflammation and innate immunity",

  "C6_Intercellular_communication",
  "Intercellular communication and glial niche",

  "C7_Retinal_cell_function",
  "Retinal cell identity and specialized function",

  "C8_Regeneration_plasticity",
  "Regeneration and developmental plasticity",

  "C9_ECM_vascular_architecture",
  "ECM, vascular function and tissue architecture"
)


## ============================================================
## 2. Curated Human Hallmark gene sets
## ============================================================

hallmark_definition <- tribble(
  ~gene_set, ~category, ~subcategory,

  ## C1: Genome and epigenome
  "HALLMARK_DNA_REPAIR",
  "C1_Genome_epigenome",
  "DNA repair",

  "HALLMARK_UV_RESPONSE_UP",
  "C1_Genome_epigenome",
  "DNA damage response",

  "HALLMARK_UV_RESPONSE_DN",
  "C1_Genome_epigenome",
  "DNA damage response",

  ## C2: Proteostasis and autophagy
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "C2_Proteostasis_autophagy",
  "Unfolded protein response",

  "HALLMARK_PROTEIN_SECRETION",
  "C2_Proteostasis_autophagy",
  "Protein processing and secretion",

  ## C3: Metabolism and mitochondria
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "C3_Metabolism_mitochondria",
  "Oxidative phosphorylation",

  "HALLMARK_GLYCOLYSIS",
  "C3_Metabolism_mitochondria",
  "Glycolysis",

  "HALLMARK_FATTY_ACID_METABOLISM",
  "C3_Metabolism_mitochondria",
  "Fatty acid metabolism",

  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "C3_Metabolism_mitochondria",
  "Cholesterol homeostasis",

  "HALLMARK_PEROXISOME",
  "C3_Metabolism_mitochondria",
  "Peroxisomal metabolism",

  "HALLMARK_MTORC1_SIGNALING",
  "C3_Metabolism_mitochondria",
  "Nutrient sensing",

  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "C3_Metabolism_mitochondria",
  "PI3K-AKT-mTOR signaling",

  "HALLMARK_XENOBIOTIC_METABOLISM",
  "C3_Metabolism_mitochondria",
  "Detoxification metabolism",

  "HALLMARK_HEME_METABOLISM",
  "C3_Metabolism_mitochondria",
  "Heme metabolism",

  ## C4: Stress, senescence and cell death
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "C4_Stress_senescence_cell_death",
  "Oxidative stress",

  "HALLMARK_P53_PATHWAY",
  "C4_Stress_senescence_cell_death",
  "P53 stress response",

  "HALLMARK_APOPTOSIS",
  "C4_Stress_senescence_cell_death",
  "Apoptosis",

  "HALLMARK_HYPOXIA",
  "C4_Stress_senescence_cell_death",
  "Hypoxia",

  ## C5: Inflammation and immunity
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "C5_Inflammation_immunity",
  "Inflammatory response",

  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "C5_Inflammation_immunity",
  "TNF-NF-kB signaling",

  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "C5_Inflammation_immunity",
  "Type I interferon response",

  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "C5_Inflammation_immunity",
  "Interferon-gamma response",

  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "C5_Inflammation_immunity",
  "IL6-JAK-STAT3 signaling",

  "HALLMARK_COMPLEMENT",
  "C5_Inflammation_immunity",
  "Complement",

  "HALLMARK_ALLOGRAFT_REJECTION",
  "C5_Inflammation_immunity",
  "Immune activation",

  "HALLMARK_IL2_STAT5_SIGNALING",
  "C5_Inflammation_immunity",
  "Cytokine signaling",

  ## C6: Intercellular communication
  "HALLMARK_NOTCH_SIGNALING",
  "C6_Intercellular_communication",
  "Notch signaling",

  "HALLMARK_TGF_BETA_SIGNALING",
  "C6_Intercellular_communication",
  "TGF-beta signaling",

  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "C6_Intercellular_communication",
  "Wnt-beta-catenin signaling",

  "HALLMARK_HEDGEHOG_SIGNALING",
  "C6_Intercellular_communication",
  "Hedgehog signaling",

  ## C8: Regeneration and plasticity
  "HALLMARK_E2F_TARGETS",
  "C8_Regeneration_plasticity",
  "Cell-cycle activation",

  "HALLMARK_G2M_CHECKPOINT",
  "C8_Regeneration_plasticity",
  "Cell-cycle activation",

  "HALLMARK_MYC_TARGETS_V1",
  "C8_Regeneration_plasticity",
  "Growth and biosynthesis",

  "HALLMARK_MYC_TARGETS_V2",
  "C8_Regeneration_plasticity",
  "Growth and biosynthesis",

  "HALLMARK_MITOTIC_SPINDLE",
  "C8_Regeneration_plasticity",
  "Cell division",

  ## C9: ECM, vascular function and architecture
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "C9_ECM_vascular_architecture",
  "ECM remodeling",

  "HALLMARK_APICAL_JUNCTION",
  "C9_ECM_vascular_architecture",
  "Cell junction",

  "HALLMARK_ANGIOGENESIS",
  "C9_ECM_vascular_architecture",
  "Angiogenesis",

  "HALLMARK_COAGULATION",
  "C9_ECM_vascular_architecture",
  "Vascular and extracellular environment"
)


## ============================================================
## 3. Curated Human GO Biological Process terms
## ============================================================

## ============================================================
## Human GO Biological Process sets: compact version
## ============================================================

go_terms <- list(

  C1_Genome_epigenome = c(
    "DNA repair",
    "cellular response to DNA damage stimulus",
    "chromatin organization",
    "histone modification",
    "nucleosome organization",
    "regulation of gene expression, epigenetic",
    "RNA processing",
    "RNA splicing",
    "telomere maintenance"
  ),

  C2_Proteostasis_autophagy = c(
    "autophagy",
    "macroautophagy",
    "protein folding",
    "response to unfolded protein",
    "endoplasmic reticulum unfolded protein response",
    "lysosome organization",
    "lysosomal transport",
    "proteasomal protein catabolic process",
    "ubiquitin-dependent protein catabolic process"
  ),

  C3_Metabolism_mitochondria = c(
    "oxidative phosphorylation",
    "cellular respiration",
    "mitochondrial respiratory chain complex assembly",
    "ATP metabolic process",
    "glucose metabolic process",
    "fatty acid metabolic process",
    "lipid metabolic process",
    "cholesterol metabolic process",
    "NAD metabolic process",
    "response to nutrient levels",
    "cellular redox homeostasis"
  ),

  C4_Stress_senescence_cell_death = c(
    "response to oxidative stress",
    "response to reactive oxygen species",
    "cellular senescence",
    "apoptotic process",
    "intrinsic apoptotic signaling pathway",
    "programmed cell death",
    "neuron apoptotic process",
    "response to hypoxia"
  ),

  C5_Inflammation_immunity = c(
    "inflammatory response",
    "innate immune response",
    "cytokine-mediated signaling pathway",
    "response to interferon-gamma",
    "type I interferon signaling pathway",
    "complement activation",
    "antigen processing and presentation",
    "phagocytosis",
    "leukocyte activation"
  ),

  C6_Intercellular_communication = c(
    "cell-cell signaling",
    "response to growth factor",
    "fibroblast growth factor receptor signaling pathway",
    "BMP signaling pathway",
    "Notch signaling pathway",
    "Wnt signaling pathway",
    "ephrin receptor signaling pathway",
    "transforming growth factor beta receptor signaling pathway"
  ),

  C7_Retinal_cell_function = c(
    "visual perception",
    "phototransduction",
    "sensory perception of light stimulus",
    "photoreceptor cell maintenance",
    "photoreceptor outer segment organization",
    "cilium organization",
    "synaptic signaling",
    "chemical synaptic transmission",
    "axonogenesis",
    "neuron projection organization",
    "neurotransmitter transport",
    "glutamate metabolic process",
    "ion homeostasis",
    "retinoid metabolic process"
  ),

  C8_Regeneration_plasticity = c(
    "tissue regeneration",
    "neurogenesis",
    "regulation of neurogenesis",
    "stem cell differentiation",
    "neural precursor cell proliferation",
    "cell cycle",
    "DNA replication",
    "cell population proliferation",
    "cell fate commitment",
    "wound healing",
    "developmental growth"
  ),

  C9_ECM_vascular_architecture = c(
    "extracellular matrix organization",
    "cell adhesion",
    "cell junction organization",
    "actin cytoskeleton organization",
    "basement membrane organization",
    "angiogenesis",
    "vasculature development",
    "blood vessel remodeling",
    "collagen metabolic process"
  )
)

go_definition <- purrr::imap_dfr(
  go_terms,
  ~ tibble::tibble(
    category = .y,
    subcategory = .x,
    go_term = .x
  )
)



########### 5. Build the unified Human aging gene-set catalog
########### 
db_species = "HS"
species = "Homo sapiens"
############

## Convert a GO term name into the MSigDB GOBP name
to_gobp_name <- function(x) {
  paste0(
    "GOBP_",
    toupper(gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", x)))
  )
}

build_human_aging_catalog <- function(hallmark_definition, go_definition) {

  ## 1. Human Hallmark gene sets
  h_all <- msigdbr::msigdbr(
    db_species = "HS",
    species = "Homo sapiens",
    collection = "H"
  )

  h <- h_all |>
    dplyr::inner_join(
      hallmark_definition,
      by = c("gs_name" = "gene_set")
    ) |>
    dplyr::transmute(
      term_id = paste0("HALLMARK::", gs_name),
      source = "Hallmark",
      gene_set = gs_name,
      category,
      subcategory,
      gene_symbol,
      db_version
    )

  ## 2. Human GO Biological Process
  g_all <- msigdbr::msigdbr(
    db_species = "HS",
    species = "Homo sapiens",
    collection = "C5",
    subcollection = "GO:BP"
  ) |>
    dplyr::mutate(
      go_id = stringr::str_extract(
        paste(
          dplyr::coalesce(gs_exact_source, ""),
          dplyr::coalesce(gs_url, ""),
          dplyr::coalesce(gs_description, "")
        ),
        "GO:[0-9]{7}"
      )
    )

  ## 3. Map the input GO terms to GO IDs
  go_key <- AnnotationDbi::select(
    GO.db,
    keys = AnnotationDbi::keys(GO.db, keytype = "GOID"),
    columns = c("TERM", "ONTOLOGY"),
    keytype = "GOID"
  ) |>
    tibble::as_tibble() |>
    dplyr::filter(ONTOLOGY == "BP") |>
    dplyr::mutate(term_key = tolower(TERM)) |>
    dplyr::inner_join(
      go_definition |>
        dplyr::mutate(term_key = tolower(go_term)),
      by = "term_key"
    ) |>
    dplyr::transmute(
      go_id = GOID,
      category,
      subcategory,
      go_term
    ) |>
    dplyr::distinct()

  g <- g_all |>
    dplyr::inner_join(go_key, by = "go_id") |>
    dplyr::transmute(
      term_id = paste0("GO_BP::", go_id),
      source = "GO_BP",
      gene_set = gs_name,
      category,
      subcategory,
      gene_symbol,
      db_version
    )

  ## 4. Check for gene sets that failed to match
  missing_hallmark <- hallmark_definition |>
    dplyr::anti_join(
      h_all |> dplyr::distinct(gs_name),
      by = c("gene_set" = "gs_name")
    )

  missing_go <- go_definition |>
    dplyr::mutate(term_key = tolower(go_term)) |>
    dplyr::anti_join(
      go_key |>
        dplyr::mutate(term_key = tolower(go_term)) |>
        dplyr::distinct(term_key),
      by = "term_key"
    )

  if (nrow(missing_hallmark) > 0) {
    warning("The following Hallmark gene sets were not matched:")
    print(missing_hallmark)
  }

  if (nrow(missing_go) > 0) {
    warning("The following GO terms were not matched in GO.db:")
    print(missing_go)
  }

  ## 5. Merge
  long <- dplyr::bind_rows(h, g) |>
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") |>
    dplyr::distinct()

  annotation <- long |>
    dplyr::count(
      term_id,
      source,
      gene_set,
      category,
      subcategory,
      db_version,
      name = "n_genes"
    )

  term2gene <- long |>
    dplyr::distinct(term_id, gene_symbol)

  ## Note: split(.$category) cannot be used here
  pathways <- split(
    term2gene$gene_symbol,
    term2gene$term_id
  ) |>
    lapply(unique)

  category_unions <- split(
    long$gene_symbol,
    long$category
  ) |>
    lapply(unique)

  list(
    annotation = annotation,
    term2gene = term2gene,
    pathways = pathways,
    category_unions = category_unions,
    missing_hallmark = missing_hallmark,
    missing_go = missing_go,
    db_version = unique(long$db_version)
  )
}

######
###### Let's inspect this ###########
######

human_aging_catalog <- build_human_aging_catalog(
  hallmark_definition = hallmark_definition,
  go_definition = go_definition
)

#######
#######

## C1: remove the obsolete "histone modification",
## replace it with "chromatin remodeling"
go_terms$C1_Genome_epigenome <- go_terms$C1_Genome_epigenome |>
  setdiff(c(
    "histone modification",
    "regulation of gene expression, epigenetic"
  )) |>
  c(
    "chromatin remodeling",
    "epigenetic regulation of gene expression"
  )

## C7: use the official GO term name
go_terms$C7_Retinal_cell_function <- go_terms$C7_Retinal_cell_function |>
  setdiff("photoreceptor outer segment organization") |>
  c("photoreceptor cell outer segment organization")


####### Save the catalog #########
#######

setwd("/projects/hmz-aging/HMZ_GO_terms_Final")

save(human_aging_catalog,file="human_aging_catalog")
########
########
######## Load #####
########
setwd("/projects/hmz-aging/HMZ_GO_terms_Final")
load("Reactome_senescence_list")
load("SenMayo_senescence_list")


senmayo_human = SenMayo_senescence_list$H
reactome_cellular_senescence_human = Reactome_senescence_list$H


## 1. Assemble the two existing Human gene sets
c10_gene_sets <- list(
  "CUSTOM::SENMAYO" = senmayo_human,
  "CUSTOM::REACTOME_CELLULAR_SENESCENCE" =
    reactome_cellular_senescence_human
) |>
  lapply(function(x) {
    unique(as.character(x[!is.na(x) & x != ""]))
  })

  


setwd("/projects/hmz-aging/HMZ_GO_terms_Final")
save(human_aging_catalog,file="human_aging_catalog")


######
###### Export each category and its subcategories to an Excel file
######

library(openxlsx)

## Category order taken directly from the loaded catalog (includes C10),
## sorted by the numeric part so C10 comes after C9 (not after C1)
all_cats <- names(human_aging_catalog$category_unions)
cat_num  <- as.numeric(sub("^C([0-9]+)_.*", "\\1", all_cats))
aging_category_order <- all_cats[order(cat_num)]

## Human-readable labels (大类); covers C1-C10
aging_category_labels <- tibble::tribble(
  ~category, ~category_label,
  "C1_Genome_epigenome",             "Genome and epigenome maintenance",
  "C2_Proteostasis_autophagy",       "Proteostasis, autophagy and lysosome",
  "C3_Metabolism_mitochondria",      "Mitochondrial metabolism and nutrient sensing",
  "C4_Stress_senescence_cell_death", "Cellular stress, senescence and cell death",
  "C5_Inflammation_immunity",        "Inflammation and innate immunity",
  "C6_Intercellular_communication",  "Intercellular communication and glial niche",
  "C7_Retinal_cell_function",        "Retinal cell identity and specialized function",
  "C8_Regeneration_plasticity",      "Regeneration and developmental plasticity",
  "C9_ECM_vascular_architecture",    "ECM, vascular function and tissue architecture",
  "C10_Cellular_senescence",         "Curated cellular senescence signatures"
)

## Subcategories (小类) from the catalog annotation
ann_tbl <- human_aging_catalog$annotation |>
  dplyr::distinct(category, subcategory)

## Fallback: any category present in the catalog but missing from annotation
## (e.g. C10 if it was only added to category_unions/pathways)
missing_cats <- setdiff(aging_category_order, unique(ann_tbl$category))
if (length(missing_cats) > 0) {
  ## known subcategory labels for categories not described in annotation
  fallback_subcats <- list(
    C10_Cellular_senescence = c("SenMayo", "Reactome cellular senescence")
  )
  fb <- purrr::map_dfr(missing_cats, function(cc) {
    subs <- fallback_subcats[[cc]]
    if (is.null(subs)) subs <- NA_character_
    tibble::tibble(category = cc, subcategory = subs)
  })
  ann_tbl <- dplyr::bind_rows(ann_tbl, fb)
}

## Build a tidy table: category label (大类) + subcategory (小类)
category_subcategory_tbl <- ann_tbl |>
  dplyr::left_join(aging_category_labels, by = "category") |>
  dplyr::mutate(
    category = factor(category, levels = aging_category_order)
  ) |>
  dplyr::arrange(category, subcategory) |>
  dplyr::transmute(
    category_id    = as.character(category),
    category_label = category_label,   # 大类
    subcategory    = subcategory        # 小类
  )

out_xlsx <- "/projects/hmz-aging/HMZ_GO_terms_Final/Aging_category_subcategory.xlsx"
openxlsx::write.xlsx(category_subcategory_tbl, file = out_xlsx, overwrite = TRUE)
message("Saved: ", out_xlsx, " (", nrow(category_subcategory_tbl), " rows)")






