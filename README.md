## Iron acquisition systems in cheese microbial communities
This repository contains scripts and Hidden Markov Models (HMMs) to detect iron acquisition genes in bacteria and fungi. These scripts and HMMs were used in this article :

Tabuteau S, Hervé V, Irlinger F, Monnet C. Metagenomic profiling and genome-centric analysis reveal iron acquisition systems in cheese-associated bacteria and fungi
DOI: [DOI]

The study investigates the iron acquisition systems (direct iron import and siderophore synthesis and import) in cheese microbial communities using metagenomes, metagenome-assembled genomes (MAGs), and genomes from the MetaPDOCheese project (Gardon et al., 2025). All details are presented in these two publications.

### Homology reserch of iron acquisition genes & HMMs design

To identify genes involved in iron acquisition systems, 148 HMMs were collected from the KOfam database (Aramaki et al., 2020) (version 2022-03-01), FeGenie (Garber et al., 2020), and Protein Family Models from the NCBI database (Li et al., 2021).  
35 HMMs were designed to enhance this set of HMMs. These HMMs are available in `HMMs/` and the script used to design them is available in `code/`.  
The set of HMMs was queried against the predicted CDSs of metagenomes, MAGs, and genomes using HMMer (https://github.com/EddyRivasLab/hmmer).  
The different scripts to run HMMer and process the data on metagenomes and MAGs/genomes are available in `code/`.

### Fungal MAGs

Fungal MAGs were retrieved from non-prokaryotic bins of the MetaPDOChese metagenomic dataset through identification of fungal markers using BUSCO (Manni et al., 2021). Taxonomic assignment of these MAGs was then performed using a phylogenetic approach with UFCG (Kim et al., 2023). The script for fungal identification and analysis is available in `code/`.

## Repository contents

- `code/` – scripts used to detect iron acquisition genes and to create HMMs in this study.
- `HMMs/` – HMM profiles to detect iron acquisition genes. Either retrieved from the literature or newly built for this project.

## Citation

If you use these scripts or HMMs, please cite:

Tabuteau S, Hervé V, Irlinger F, Monnet C. Metagenomic profiling and genome-centric analysis reveal iron acquisition systems in cheese-associated bacteria and fungi
DOI: [DOI]

## References:

(Gardon et al., 2025)

Aramaki, T., Blanc-Mathieu, R., Endo, H., Ohkubo, K., Kanehisa, M., Goto, S., and Ogata, H. (2020) KofamKOALA: KEGG Ortholog assignment based on profile HMM and adaptive score threshold. Bioinformatics 36: 2251–2252.

Garber, A.I., Nealson, K.H., Okamoto, A., McAllister, S.M., Chan, C.S., Barco, R.A., and Merino, N. (2020) FeGenie: A Comprehensive Tool for the Identification of Iron Genes and Iron Gene Neighborhoods in Genome and Metagenome Assemblies. Front Microbiol 11: 37.

Li, W., O’Neill, K.R., Haft, D.H., DiCuccio, M., Chetvernin, V., Badretdin, A., et al. (2021) RefSeq: expanding the Prokaryotic Genome Annotation Pipeline reach with protein family model curation. Nucleic Acids Research 49: D1020–D1028.

Manni, M., Berkeley, M.R., Seppey, M., Simão, F.A., and Zdobnov, E.M. (2021) BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. Molecular Biology and Evolution 38: 4647–4654.

Kim, D., Gilchrist, C.L.M., Chun, J., and Steinegger, M. (2023) UFCG: database of universal fungal core genes and pipeline for genome-wide phylogenetic analysis of fungi. Nucleic Acids Research 51: D777–D784.

## License
This project is licensed under the GNU General Public License V3 - see the LICENSE file for details.
