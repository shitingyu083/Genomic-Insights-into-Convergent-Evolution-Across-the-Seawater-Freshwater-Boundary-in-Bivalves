# Policy and Method Basis for the Code Package

This file records the external basis used to audit the workflow structure.

## Nature and FAIR Requirements

- Nature Portfolio requires authors to describe how custom code and software supporting the conclusions are available, and to provide sufficient information for reuse and reproducibility where possible.
- Springer Nature research-data policy emphasizes availability of data supporting the findings, clear repository or accession information, and source data for figures where applicable.
- FAIR principles require that data and metadata be findable, accessible, interoperable and reusable.
- Software-citation principles recommend clear identification of software, versioning, authorship and persistent access.

How the workflow responds:

- `README.md` explains scope, inputs, dependencies and interpretation boundaries.
- `METHOD_CODE_MAP.csv` links each manuscript component to scripts, inputs and outputs.
- `external/external_software_commands.md` records upstream command templates.
- `data/templates/` provides machine-readable input schemas.
- `CODE_FILE_MANIFEST.csv` records package contents.
- `outputs/sessionInfo.txt` is written at runtime.

## Statistical and Bioinformatic Method Basis

- PAML/codeml branch-site outputs are summarized with both conservative chi-square(1) and 50:50 mixture-null P values for the modified branch-site framework.
- PGLS analyses use Pagel's lambda correlation structure through `ape` and `nlme`, matching the manuscript's phylogenetically informed comparative modelling.
- Gene-family copy-number screening uses `phytools::phylANOVA` for phylogenetic ANOVA.
- Amino-acid composition is treated as compositional data through centred log-ratio transformation, followed by phylogenetically informed tests and BH-FDR correction.
- TE Kimura divergence is handled as a relative divergence-from-consensus proxy, not an absolute insertion-time estimate.

## Key References to Cite in Documentation or Methods

- Wilkinson, M. D. et al. The FAIR Guiding Principles for scientific data management and stewardship. Sci. Data 3, 160018 (2016).
- Smith, A. M. et al. Software citation principles. PeerJ Comput. Sci. 2, e86 (2016).
- Yang, Z. PAML 4: phylogenetic analysis by maximum likelihood. Mol. Biol. Evol. 24, 1586-1591 (2007).
- Zhang, J., Nielsen, R. & Yang, Z. Evaluation of an improved branch-site likelihood method for detecting positive selection at the molecular level. Mol. Biol. Evol. 22, 2472-2479 (2005).
- Aitchison, J. The Statistical Analysis of Compositional Data. Chapman and Hall (1986).
- Revell, L. J. phytools: an R package for phylogenetic comparative biology and other things. Methods Ecol. Evol. 3, 217-223 (2012).
