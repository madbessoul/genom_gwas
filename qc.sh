# PROJET GWAS
# Santy
# Bessoul

# Script pour effectuer le quality control et
# l'analyse de stratification et générer des
# des fichiers de données propres

for f in $1/*.map
do

  # Récupérer le nom de base du fichier à traiter
  describer=`basename $f .map`

  # Créer un dossier pour stocker les fichiers propres
  # pour  chaque chromosome
  mkdir $1/${describer}


  # Executer la commande plink pour effectuer le quality control
  plink --noweb --file $1/${describer} --hwe 0.001 --geno 0.02 --maf 0.01 --recode --out $1/${describer}/${describer}_RC

  cd $1/${describer}

  # Faire le clustering et créer les matrices de distance
  plink --noweb --file ${describer}_RC --cluster --distance-matrix --out ${describer}_clust

  # Ajouter les identifiants de colonnes et de lignes
  # dans la matrice de distance
  Rscript ../../src/dist_mat.R ${describer}_clust.cluster2 ${describer}_clust.mdist

done
