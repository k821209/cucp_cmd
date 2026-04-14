# Castilleja Closest Phylogenetic Neighbors

## Overview

*Castilleja* (Indian paintbrush) belongs to the family **Orobanchaceae** (broomrape family), a diverse group of parasitic and hemiparasitic plants. This document lists genera in order of phylogenetic relationship to *Castilleja* for reference dataset construction.

## Phylogenetic Ranking (Closest First)

| Rank | Genus | Tribe/Clade | Relationship to Castilleja | Notes |
|------|-------|-------------|---------------------------|-------|
| 1 | ***Castilleja*** | Castillejinae | **Focal taxon** | Hemiparasitic |
| 2 | *Pedicularis* | Pedicularideae | Sister to Castillejinae | Louseworts; ~600 species |
| 3 | *Euphrasia* | Rhinantheae | Same subtribe cluster | Eyebrights |
| 4 | *Bartsia* | Rhinantheae | Same subtribe cluster | Alpine bartsia |
| 5 | *Rhinanthus* | Rhinantheae | Same subtribe cluster | Yellow rattles |
| 6 | *Lathraea* | Rhinantheae | Holoparasite in tribe | Toothworts |
| 7 | *Melampyrum* | Rhinantheae | Basal Rhinantheae | Cow-wheats |
| 8 | *Triaenophora* | Basal Orobanchaceae | Basal hemiparasite | Monotypic genus |
| 9 | *Rehmannia* | Basal Orobanchaceae | Basal clade | Chinese foxglove |
| 10 | *Lindenbergia* | Basal Orobanchaceae | Non-parasitic basal | Sister to parasitic clade |

## Outgroups

| Rank | Genus | Family | Relationship |
|------|-------|--------|--------------|
| 11 | *Erythranthe* | Phrymaceae | Close outgroup (Lamiales) |
| 12 | *Mimulus* | Phrymaceae | Close outgroup (Lamiales) |
| 13 | *Paulownia* | Paulowniaceae | Distant outgroup (Lamiales) |

## Excluded Genera

### Holoparasites with Reduced Plastomes

These genera were excluded due to significant plastome reduction (<140 kb) causing alignment artifacts and long-branch attraction:

| Genus | Plastome Size | Reason |
|-------|--------------|--------|
| *Orobanche* | ~88 kb | Holoparasite - reduced plastome |
| *Phelipanche* | ~63 kb | Holoparasite - reduced plastome |
| *Cistanche* | ~99 kb | Holoparasite - reduced plastome |
| *Epifagus* | ~70 kb | Holoparasite - reduced plastome |
| *Conopholis* | ~46 kb | Holoparasite - reduced plastome |
| *Boulardia* | ~80 kb | Holoparasite - reduced plastome |

### Structural Anomalies

| Genus | Plastome Size | Reason |
|-------|--------------|--------|
| *Striga* | ~191 kb | Massive IR expansion (~70 kb vs typical ~25 kb) |

## Sequence Counts (in Reference Dataset)

| Genus | Representative | Max 3 per genus |
|-------|----------------|-----------------|
| *Castilleja* | 1 | 1 |
| *Pedicularis* | 45 | 3 |
| *Euphrasia* | 2 | 2 |
| *Bartsia* | 1 | 1 |
| *Rhinanthus* | 3 | 3 |
| *Lathraea* | 2 | 2 |
| *Melampyrum* | 4 | 3 |
| *Triaenophora* | 1 | 1 |
| *Rehmannia* | 7 | 3 |
| *Lindenbergia* | 2 | 2 |
| *Erythranthe* | 18 | 3 |
| *Mimulus* | 2 | 2 |
| *Paulownia* | 9 | 3 |
| **TOTAL** | **97** | **29** |

## Key References

1. **Tank & Olmstead (2008)** "From annuals to perennials: phylogeny of subtribe Castillejinae (Orobanchaceae)" *American Journal of Botany*

2. **McNeal et al. (2013)** "Phylogeny and origins of holoparasitism in Orobanchaceae" *American Journal of Botany*

3. **Yu et al. (2018)** "Phylogenetic Relationships in Orobanchaceae Inferred From Low-Copy Nuclear Genes" *Frontiers in Plant Science* https://pmc.ncbi.nlm.nih.gov/articles/PMC6646720/

4. **Frailey et al. (2018)** "Gene loss and genome rearrangement in the plastids of five Hemiparasites in the family Orobanchaceae" *BMC Plant Biology* 18:30

## Reference Files

| File | Sequences | Description |
|------|-----------|-------------|
| `Orobanchaceae_cp_genomes_representative.fasta` | 97 | One sequence per species |
| `Orobanchaceae_cp_genomes_max3.fasta` | 29 | Max 3 species per genus |
| `Castilleja_cp_complete_genomes.fasta` | 2 | Only *Castilleja* sequences |
