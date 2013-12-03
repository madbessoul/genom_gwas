#! /bin/bash

# PROJET GWAS
# Santy
# Bessoul

# Script pour effectuer le quality control et
# l'analyse de stratification et générer des
# des fichiers de données propres

# Le script traite automatiqement tous les chromosomes
# dans un dossier donné (donné en premier paramètres)


for f in $1/*.map
do

  # Récupérer le nom de base du fichier à traiter
  describer=`basename $f .map`

  # Créer un dossier pour stocker les résultats
  # pour  chaque chromosome
  mkdir $1/${describer}

  echo 'Quality control'
  echo '==============='
  # Executer la commande plink pour effectuer le quality control
  plink --noweb --file $1/${describer} --hwe 0.001 --geno 0.02 --maf 0.01 --recode --out $1/${describer}/${describer}_RC

  echo Stratification Analysis
  echo =======================
  # Faire le clustering et créer les matrices de distance
  plink --noweb --file $1/${describer}/${describer}_RC --cluster --distance-matrix --out $1/${describer}/${describer}_clust

  # Retrait des outliers et sauvegarde des individus à enlever
  # dans toremove.txt

  Rscript dist_mat.R $1/${describer}/${describer}_clust.cluster2 $1/${describer}/${describer}_clust.mdist

  mv toremove.txt $1/${describer}
  mv *.pdf $1/${describer}

  # Génerer des fichiers propres en enlever les outliers des data
  plink --noweb --file $1/${describer}/${describer}_RC --remove $1/${describer}/toremove.txt --recode --out $1/${describer}/${describer}_QC

  echo Removing intermediate data files
  echo ================================
  rm -rf $1/${describer}/${describer}_RC.*


done
