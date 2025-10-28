---
title: "R Notebook"
output: github_document
---
```{r}
```

 Le chargement du package dada2
# Ce package est utilisé pour analyser les séquences d’ADN issues du séquençage haut débit (NGS),
# notamment celles du gène 16S dans ce travail. Il permet de filtrer, corriger, assembler
# et classifier les séquences afin d’obtenir des ASVs (Amplicon Sequence Variants)
# avec une meilleure précision que les méthodes traditionnelles.
# Chargement du package dada2

```{r}
library(dada2)
```

# Vérification de la version du package dada2 installée
packageVersion("dada2")

```{r}
library(dada2); packageVersion("dada2")
```

# Définition du chemin où se trouvent mes fichiers fastq (les données brutes de séquençage Illumina)
# Ici, le chemin est "~/MiSeq_SOP" : 
#   - "~" correspond à mon répertoire utilisateur (home directory).
#   - "MiSeq_SOP" est le dossier qui contient mes fichiers .fastq après décompression.

```{r}
setwd("/home/rstudio/tuturiol_sow")
path <- "~/MiSeq_SOP" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
```


# Les fichiers FASTQ forward et reverse suivent un format de nommage standard :
# Exemple : SAMPLENAME_R1_001.fastq et SAMPLENAME_R2_001.fastq
#   - R1 = lectures forward (lecture avant)
#   - R2 = lectures reverse (lecture arrière)
#   - "SAMPLENAME" = identifiant de ton échantillon

# On crée une liste des fichiers forward (R1)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))

# On crée une liste des fichiers reverse (R2)
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extraction des noms d'échantillons
# basename(fnFs) → garde seulement le nom du fichier (pas le chemin complet)
# strsplit(..., "_") → découpe le nom du fichier en morceaux séparés par "_"
# sapply(..., `[`, 1) → prend le 1er morceau (donc le nom de l'échantillon avant le "_")
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

```{r}
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```


# Visualisation des profils de qualité des lectures Forward
# Ici, on affiche les graphiques de qualité pour les deux premiers fichiers FASTQ "Forward".
# Cette analyse permet de comparer la dégradation de la qualité en fin de lecture
# et d’ajuster les paramètres de trimming en conséquence.

```{r}
plotQualityProfile(fnFs[1:2])
```



# Visualisation des profils de qualité des lectures Reverse
# Comme pour les reads Forward, on inspecte ici la qualité des deux premiers fichiers FASTQ "Reverse".
# Cette analyse permet de comparer la dégradation de la qualité en fin de lecture
# et d’ajuster les paramètres de trimming en conséquence.

```{r}
plotQualityProfile(fnRs[1:2])
```
# Création des chemins de sortie pour les fichiers filtrés (lectures forward)
# On utilise file.path() pour faire un chemin complet :
#   - "path" = dossier principal (~/MiSeq_SOP)
#   - "filtered" = sous-dossier où seront stockés l es fichiers filtrés
#   - paste0() le nom de l’échantillon avec le suffixe "_F_filt.fastq.gz"
# meme chose pour les lectures reverse (R2), suffixe "_R_filt.fastq.gz"
# On associe chaque fichier forward filtré à son nom d’échantillon
names(filtFs) <- sample.names
# Idem pour les fichiers reverse filtrés
names(filtRs) <- sample.names

```{r}
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

# Filtrage et tronquage des lectures forward et reverse
out <- filterAndTrim(
#fnFs,fichiers FASTQ forward d’origine
#filtFs,chemins de sortie pour forward filtrés
#fnRs,  fichiers FASTQ reverse d’origine
# filtRs, chemins de sortie pour reverse filtrés
# truncLen = c(240,160),  Longueurs de tronquage : 240 bases pour R1, 160 pour R2
# maxN = 0, Aucun N autorisé (bases indéterminées)
#maxEE = c(2,2),Maximum attendu d’erreurs : 2 pour forward, 2 pour reverse
#truncQ = 2, Tronquer à partir d’une base < Q2 (très basse qualité)
#rm.phix = TRUE,Supprime les lectures correspondant au contrôle PhiX
#compress = TRUE, Fichiers de sortie compressés en .gz
#multithread = F, Sur Windows mettre FALSE, sur Linux/Mac on peut mettre TRUE

```{r}
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=F) # On Windows set multithread=FALSE
head(out)
```
# Apprentissage du modèle d’erreurs pour les lectures Forward filtrées
# Cette étape permet à dada2 d’estimer les profils d’erreurs propres aux données,
# afin de distinguer les erreurs de séquençage des variations biologiques réelles.

```{r}
errF <- learnErrors(filtFs, multithread=TRUE)
```

# Apprentissage du modèle d’erreurs pour les lectures Reverse filtrées
# Comme pour les reads Forward, dada2 estime ici les erreurs spécifiques aux données Reverse,
# ce qui permettra une correction plus fiable des séquences lors de l’inférence des ASVs.

```{r}
errR <- learnErrors(filtRs, multithread=TRUE)
```


# Visualisation du modèle d’erreurs des lectures Forward
# Cette représentation permet de vérifier si le modèle d’erreur estimé
# correspond bien aux tendances attendues (en fonction des scores de qualité).

```{r}
plotErrors(errF, nominalQ=TRUE)
```




# Inférence des variants de séquences exactes (ASVs) à partir des lectures Forward filtrées
# Cette étape applique le modèle d’erreurs pour corriger les séquences et identifier
# les variants uniques réels. Le résultat indique notamment le nombre de lectures conservées
# et le nombre d’ASVs détectés.

```{r}
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```


# Inférence des variants de séquences exactes (ASVs) à partir des lectures Reverse filtrées
# Comme pour les séquences Forward, le modèle d’erreurs est utilisé ici pour corriger
# les lectures Reverse et identifier les ASVs réellement présents dans l’échantillon.

```{r}
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```



# Consultation du résultat pour le premier échantillon
# Utile pour obtenir un aperçu du nombre de séquences uniques (ASVs “bruts”)
# identifiées dans cet échantillon avant la fusion des lectures.
```{r}
dadaFs[[1]]
```


# Fusion des lectures Forward et Reverse (R1 et R2)
# Cette étape assemble les paires de lectures qui se chevauchent afin de reconstruire
# la séquence complète de l’amplicon. On affiche ici les premières lignes du résultat
# correspondant au premier échantillon pour vérifier la qualité de la fusion.

```{r}
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```



# Construction de la table des variants de séquences (ASV table)
# Cette matrice contient, pour chaque échantillon, le nombre de lectures
# correspondant à chaque ASV obtenu après fusion des lectures.
# Affichage des dimensions de la table : nombre d’échantillons x nombre d’ASVs

```{r}
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```


# Vérification de la distribution des longueurs des séquences (ASVs)
# Cette commande génère une table de fréquence indiquant combien d’ASVs correspondent à chaque longueur de séquence, Cela permet de repérer
# d’éventuelles séquences aberrantes avant la suite du traitement.


```{r}
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
```


# Suppression des séquences chimériques
# Dimensions de la nouvelle table (après suppression des chimères)
#Détecte les chimères = artefacts créés lors de la PCR, où deux séquences différentes s’assemblent faussement
#C’est une étape indispensable, car les chimères peuvent représenter jusqu’à 20–30 % des séquences dans un jeu 16S
#Resultat, Table d’ASVs nettoyée : seules les séquences biologiquement plausibles restent

```{r}
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
```



# Calcul de la proportion de lectures non chimériques
# On compare le nombre total de lectures après suppression des chimères
# au nombre total initial, afin d’évaluer la fraction de données réellement conservées.

```{r}
sum(seqtab.nochim)/sum(seqtab)
```


# Fonction utilitaire : compte le nombre total de lectures uniques
# Construction du tableau de suivi
# Attribution des noms de colonnes
# Attribution des noms de lignes (noms d’échantillons)
# Affiche les 6 premières lignes du tableau
#Ça permet de visualiser les pertes à chaque étape (filtrage, débruitage, fusion, anti-chimères).
#On peut rapidement repérer un problème (ex : très faible taux de fusion ou énorme perte après suppression des chimères).

```{r}
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```




# Création d’un dossier “tax” dans le répertoire utilisateur (home directory)
# Ce dossier contiendra les fichiers de référence nécessaires à l’assignation taxonomique.

# Téléchargement / mise en place du jeu de référence Silva v132

## Base qui permet l’assignation taxonomique jusqu’au niveau du genre (assignTaxonomy)

# Téléchargement du fichier complémentaire Silva v132

## Utilisé pour l’assignation jusqu’au niveau espèce lorsque cela est possible (addSpecies)


```{r}
dir.create("~/tax", showWarnings=FALSE)
download.file("https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz",
              destfile = "~/tax/silva_nr_v132_train_set.fa.gz",mode = "wb")
download.file("https://zenodo.org/records/1172783/files/silva_species_assignment_v132.fa.gz",
              destfile = "~/tax/files/silva_species_assignment_v132.fa.gz",mode = "wb")
            
```



# Assignation taxonomique des ASVs à partir de la base de référence SILVA v132
# Cette étape associe à chaque ASV une classification taxonomique
# du niveau Domaine jusqu’au Genre (selon le degré de similarité des séquences).
```{r}
taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
```


# Ajout de l’assignation taxonomique au niveau espèce
# On utilise ici le fichier de référence “silva_species_assignment_v132.fa.gz”
# contenant des séquences annotées au niveau espèce, afin de compléter
# la classification taxonomique lorsque cela est possible.
```{r}
taxa <- addSpecies(taxa, "~/silva_species_assignment_v132.fa.gz")
```



# Préparation d’une version affichable de la table taxonomique
# On supprime les noms de lignes (séquences ASVs) pour faciliter la lecture
# et afficher un aperçu clair des classifications taxonomiques obtenues.
```{r}
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```



# Chargement du package DECIPHER
# Ce package est dédié à l’analyse de séquences ADN/ARN, notamment pour l’alignement
# multiple, l’identification fonctionnelle et d’autres analyses bio-informatiques avancées.
library(DECIPHER)

# Vérification de la version installée du package
packageVersion("DECIPHER")

```{r}
library(DECIPHER); packageVersion("DECIPHER")
```



# Création d’un objet contenant les séquences ASVs au format DNAStringSet
dna <- DNAStringSet(getSequences(seqtab.nochim))

# Chargement du jeu de référence DECIPHER (SILVA R138.2)
# Ce fichier .RData contient le modèle d’apprentissage nécessaire pour l’assignation taxonomique.
load("~/SILVA_SSU_r138_2_2024.RData") # À adapter selon l'emplacement du fichier

# Assignation taxonomique via DECIPHER (IdTaxa)
# Utilisation du brin “top” et de l’ensemble des processeurs disponibles
ids <- IdTaxa(dna, trainingSet, strand = "top",
              processors = NULL, verbose = TRUE)

# Définition des rangs taxonomiques d’intérêt pour la sortie finale
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

# Conversion du résultat (objet de classe "Taxa") en matrice taxonomique,
# structurée comme celle issue de assignTaxonomy()
##Les annotations “unclassified_...” sont remplacées par NA pour plus de clarté
taxid <- t(sapply(ids, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
}))

# Attribution des noms de colonnes et des ASVs comme noms de lignes
colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab.nochim)

```{r}
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
  load("~/SILVA_SSU_r138_2_2024.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=TRUE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
```



# Sélection des ASVs présents dans l’échantillon “Mock”
# On conserve uniquement les ASVs avec une abondance > 0,
# puis on les trie par abondance décroissante.

# Affichage du nombre total d’ASVs détectés dans le Mock community.
# Dans un Mock, le nombre d’espèces attendues est connu (ex. ~20 ici),
# ce qui permet d’évaluer la performance du pipeline.

```{r}
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
```

# Chargement des séquences de référence du Mock community
# Le fichier “HMP_MOCK.v35.fasta” contient les séquences 16S des espèces connues du Mock,
# servant ici de jeu de référence pour évaluer la précision de l’analyse.

# Vérification de la concordance entre les ASVs détectés et les séquences attendues
# On teste pour chaque ASV du Mock si une correspondance exacte existe dans les références.

# Affichage du nombre d’ASVs correspondant parfaitement au jeu d’espèces attendu


```{r}
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
```


# Chargement du package phyloseq
# Vérification de la version installée
#Il est conçu pour faciliter l’intégration, la manipulation et la visualisation de données issues de pipelines comme DADA2
```{r}
library(phyloseq); packageVersion("phyloseq")
```


# Chargement du package Biostrings
# Vérification de la version installée
#charge le package Biostrings, qui fournit des classes et fonctions pour manipuler efficacement des séquences biologiques (ADN, ARN, protéines).

```{r}
library(Biostrings); packageVersion("Biostrings")
```


# Chargement du package ggplot2
# Vérification de la version installée
#charge le package ggplot2, qui fait partie du tidyverse et permet de créer des graphiques élégants et personnalisés.

```{r}
library(ggplot2); packageVersion("ggplot2")
```
# Définit le thème par défaut pour tous les graphiques ggplot2
#permet de changer le thème global appliqué à tous tes graphiques ggplot2

```{r}
theme_set(theme_bw())
```



# Récupère les noms des échantillons (correspondent aux lignes de la table ASV)
# Sépare le nom de l’échantillon sur le "D" et prend la première partie (= sujet/identifiant)
# Le premier caractère du "subject" = le genre (M = male, F = female)
# Le reste du "subject" (après le 1er caractère) = identifiant du sujet
# Deuxième partie après "D" = numéro du jour (converti en entier)
# Création d’un data.frame avec les colonnes Sujet, Genre et Jour
# Ajout d’une variable "When" = précoce (Day ≤ 100) ou tardif (Day > 100)
# Les noms des lignes du data.frame sont les noms des échantillons
#ça m'a permi de créer une table de métadonnées par échantillon (sujets, genre, jour, période)

```{r}
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out
```

# Création de l’objet phyloseq regroupant :
## la table d’abondances des ASVs
##  les métadonnées associées aux échantillons
## la table taxonomique issue de l’assignation
# Suppression de l’échantillon témoin “Mock” pour les analyses écologiques ultérieures


```{r}
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
```


# Création d’un objet DNAStringSet avec les séquences des ASVs
# Attribution des noms de séquences : les ASVs sont nommés avec leurs séquences nucléotidiques complètes
# Ajout de ces séquences à l’objet phyloseq (fusion des données et des séquences)
# Renommage des ASVs avec des identifiants courts et lisibles (ASV1, ASV2, …)
# Affichage de l’objet phyloseq (résumé : nombre d’échantillons, d’ASVs, de taxons, etc.)
# Les séquences brutes étant trop longues pour servir de noms, leur remplacement par des IDs courts facilite les affichages,
# tout en conservant les séquences si besoin pour l’alignement ou la phylogénie.

```{r}
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
```


# Tracé de la diversité alpha en fonction du jour d’échantillonnage
# ps : objet phyloseq utilisé comme source de données
# x = "Day" : la variable Day est représentée en abscisse
# measures = c("Shannon", "Simpson") : indices de diversité alpha calculés et affichés
# Le graphique permet de visualiser l’évolution de la diversité microbienne au cours du temps
# et d’observer les différences selon la variable When (Early vs Late).

```{r}
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")
```

# Transformation des abondances en proportions (abondances relatives par échantillon)
# Chaque échantillon est normalisé de sorte que la somme de ses abondances soit égale à 1

# Réalisation d’une ordination NMDS basée sur la distance de Bray–Curtis
# à partir des abondances relatives

```{r}
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
```

# Trace l’ordination NMDS calculée précédemment
#plot_ordination() : fonction phyloseq pour représenter une ordination (ordinate)

```{r}
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")
```


# Identifie les 20 ASVs les plus abondants (somme sur tous les échantillons)
# Transforme les abondances en proportions par échantillon
# Garde uniquement les 20 ASVs les plus abondant
# Trace un barplot empilé par échantillon (x=Day), avec les familles en couleurs
# et facettes séparées selon la variable "When" (Early vs Late)
```{r}
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")
```



