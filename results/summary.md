# HNSC Mutation Landscape Analysis -- Results Summary


## Overview


Analyzed 528 TCGA Head and Neck Squamous Cell Carcinoma primary tumor samples, of which 510 had somatic mutation data. After excluding 6 hypermutated samples (identified by 3-standard-deviation method), the final analytical cohort comprised 504 samples with 84,690 non-silent somatic mutations. Four analyses were performed: mutation frequency ranking, oncoplot visualization, pariwise co-occurrence testing, and HPV-stratified mutation comparison. All statistical tests use Benjamini-Hochberg FDR correction (q < 0.05).


## Figure 1: Mutation Frequency

Tp53 is mutated in 70.6% of non-hypermutated samples, confirming its role as the most recurrently mutated genes in HNSC. Other frequently mutated known drivers include FAT1 (22.6%), CDKN2A (21.8%), PIK3CA (18.5%), NOTCH1 (17.9%), KMT2D (16.1%), NSD1 (11.5%), and CASP8 (10.3%). Togther, these eight genes represent the core driver landscape of HNSC, disrupting cell cycle control (TP53, CDKN2A), Wnt signaling (FAT1), PI3K signaling (PIK3CA), squamous differentiation (NOTCH1), chromatic regulation (KMT2D, NSD1), and apoptosis (CASP8).

Several large genes appear in top 24, like TTN (41.7%), CSMD3 (19.8%), MUC16 (18.8%), LRP1B (16.7%), SYNE1 (17.5%). But these are likely passengers mutated proportionally to their coding sequence length rather than by positive selection. This interpretation is supported by their mutation type profiles: predominantly missense with high silent fractions in the unfiltered data. Gene-length correction (e.g., MutSig2CV) was not performed in this analysis.

FRG1B (20%) is a probable alignment artifact arising from read mismapping among paralogous FSHD Region Gene 1 family members. It does not appear in published significantly mutated gene lists for HNSC.


## Figure 2: Oncoplot

The oncoplot displays mutation types across the top 20 genes in 482 samples (samples with at least one mutation in these genes, out of 504 non-hypermutated). Samples are sorted by TP53 mutation status.

Visually, mutation type patterns distinguish drivers from passengers. TP53 shows diverse mutation types (missense, nonsense, splice site, and frameshift) which is characterstic of a tumor suppressor undergoing loss-of-function through multiple mechanisms. FAT1 is dominated by truncating mutations (45.6% nonsense, 24.4% frameshift in the unfiltered exploration data), indicating strong selection for complete gene inactivation. PIK3CA is almost exclusively missense. TTN and other large genes show predominantly missense mutations, consistent with random accumulation.


## Figure 3: Gene Co-occurrence

Fisher's exact test was applied to all 105 pairwise combinations of the top 15 mutated genes (504 samples). After Benjamini-Hochberg correction, 15 pairs were significant at q < 0.05. All 15 showed co-occurrence; no significant mutual exclusivity was detected.

**Driver co-occurrence (3 pairs):**

TP53 and CDKN2A show the strongest co-occurrence in the dataset (OR = 12, q < 0.01). This was initially unexpected as both regulate cell cycle checkpoints, suggesting potential redundancy. However, CDKN2A encodes two distinct proteins from alternative reading frames: p16/INK4a (which inhibits CDK4/6 in the Rb pathway) and p14ARF (which stabilizes p53 via MDM2 inhibition). In a TP53-mutant tumor, p14ARF function is already irrelevant, but p16-mediated Rb pathway control remains a separate tumor-suppressive function. Thus, TP53 mutation and CDKN2A loss disable different pathways (p53 signaling and Rb signaling, respectively) explaining the co-occurrence (Sherr, 2001, Nature Reviews Molecular Cell Biology)

FAT1 and NOTCH1 co-occur (OR = 2.5, q = 0.005), consistent with tumors benefiting from disrupting both Wnt signaling (FAT1) and squamous differentiation (NOTCH1) — functionally independent pathways. FAT1 and CDKN2A also co-occur (OR = 2.1, q = 0.023).


**Passenger co-occurrence (12 pairs):**

The remaining 12 significant pairs involve combinations of large genes (TTN, CSMD3, MUC16, SYNE1, LRP1B, PCLO, FLG, DNAH5) with odds ratios between 2 and 3.3. The minimal explanation is shared mutation burden: samples with higher overall mutation counts are more likely to harbor mutations in multiple large genes simultaneously. This was not tested (e.g., by logistic regression controlling for total mutation count) and is noted as a limitation.

**TP53-PIK3CA mutual exclusivity** showed a trend (OR = 0.53, uncorrected p = 0.008) consistent with the known enrichment of PIK3CA mutations in HPV-positive HNSC, but was not observed after BH correction (q = 0.054). This is likely a power issue rather than absence of a real biological association.


## Figure 4: HPV-Stratified Mutation Comparison

Mutation frequencies were compared between HPV-positive (n = 35) and HPV-negative (n = 73) samples using Fisher's exact test with BH correction. HPV status was determined by p16 immunohistochemistry. Three genes showed significantly different mutation frequencies:

**TP53: 14.3% in HPV+ vs. HPV- (q < 0.001).** HPV oncoprotein E6 targets p53 for proteosomal degradation via the E6AP ubiquitin ligase, eliminating selective pressure for somatic TP53 mutations (Scheffner et al., 1990, Cell). The 5 HPV-positive samples with TP53 mutations may represent p16 false positives or rare HPV-positive tumors with additional somatic hits.

**FAT1: 2.9% in HPV+ vs. 26% in HPV- (q = 0.023).** Only one HPV-positive sample carries a FAT1 mutation. FAT1 loss appears characteristic of the HPV-negative, tobacco/alcohol-driven subtype.

**CDKN2A: 5.7% in HPV+ vs. 27.4% in HPV- (q = 0.049).** Borderline significant. HPV E7 binds and inactivates Rb protein directly, bypassing the need for p16 (CDKN2A) loss (Munger et al., 2004, Journal of Virology). A slightly different analytical choice could shift this result above or below the significance threshold.

**PIK3CA** trended toward enrichment in HPV-positive tumors (22.9% vs. 12.3%) consistent with published findings, but did not reach statistical significance (q = 0.51), likely due to insufficient power with 35 HPV-positive samples.


## Limitations

This analysis has several important limitations. First, mutation frequencies are not corrected for gene length or background mutation rare; large genes (TTN, CSMD3, MUC16) appear frequently mutated due to their size, not necessarily due to selection. Tools such as MutSig2CV address this but were not applied here. Second, HPV status was determined by p16 immunohistochemistry, a surrogate marker with known limitations in specificity, particularly in non-oropharyngeal sites.. Third, the HPV-positive group (n = 35 with mutation data) limits statistical power to detecting large effects only. Fourth, passenger gene co-occurrence was not formally distinguished from mutation-burden-driven correlation. Fifth, copy number alterations were not integrated -- CDKN2A is more commonly deleted than mutated in HNSC, meaning our mutation-only analysis underestimates its true inactivation frequency. Sixth, no survival or treatment response anlyses were performed; all association are descriptive.


## Future Directions

Natural extensions of this work include untegrating copy number data (available in the same cBioPortal download as GISTIC results) to capture homozygous deletions, particularly for CDKN2A. Comparison of this TCGA baseline against mutation data from South Asian HNSC patients -- where HPV prevalence and tobacco exposure patterns differ -- would address the underrepresentation of these populations in current genomic databases.