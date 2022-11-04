#!/bin/bash

########################################################################################################
# correr snippy sobre todos los ensambles en un directorio (ASSEMBLY)                                  #
# eligiendo el mejor ensamble de referencia en base a N50 y No de contigs de ser posible, es decir,    #
# que coincida el mejor ensamble en base a ambas estadisticas                                          #
# de no ser asi, se elige el mejor ensamble solo en base a No de contigs                               #
# este script requiere d las estadisticas de ensambles en /ASSEMBLY/Stats/estadisticas_ensamble.csv    #
########################################################################################################


run_snippy() {

   ##########################################
   # elegir el mejor ensamble de referencia #
   ##########################################

   # guardar valor de numero de contigs mas pequeño
   contigs_min=$(cat ASSEMBLY/Stats/estadisticas_ensamble.tsv | awk '{print $2}' | awk '/[0-9]/' | sort -n | sed -n '1p')
   # guardar valor de N50 mas grande
   N50_max=$(cat ASSEMBLY/Stats/estadisticas_ensamble.tsv | awk '{print $5}' | awk '/[0-9]/' | sort -nr | sed -n '1p')
   # obtener ensamble con numero de contigs menor
   name_contigs=$(awk -v min=${contigs_min} '$2 == min' ASSEMBLY/Stats/estadisticas_ensamble.tsv | sed -n '1p' | awk '{print $1}')
   # obtener ensamble con N50 mayor
   name_N50=$(awk -v max=${N50_max} '$5 == max' ASSEMBLY/Stats/estadisticas_ensamble.tsv | sed -n '1p' | awk '{print $1}')
   # itera sobre todos los ensambles del directorio ASSEMBLY
   for file in ASSEMBLY/*.fa; do
      # si el ensamble elegido en base a No de contigs y N50 es el mismo, entonces
      if [[ ${name_contigs} == ${name_N50} ]]; then
         # y si ademas el nombre corto del ensamble aparece en el nombre completo de ensambles de la carpeta ASSEMBLY, entonces
         if [[ ${file} =~ ${name_contigs} ]]; then
            # usar ese como ensamble de referencia
            echo "####################################################################################################################"
            echo "  el siguiente genoma se selecciono en base a los valores de N50 de: ${N50_max} y No de contigs de: ${contigs_min}  "
            echo "  se usara como genoma de referencia a = ${file}                                                                    "
            echo "####################################################################################################################"
            echo ""
            # asignar nombre clave de genoma de referencia
            ref_name=$(basename ${file} | cut -d '/' -f '2' | cut -d '.' -f '1' | cut -d '-' -f '1')
            # asignar directorio de salida para guardar resultados de Snippy
            dir="SNIPPY_${ref_name}"
         fi
      # si los nombres de ensables no son iguales, utiliza solamente el nombre del ensamble obtenido en base a No de contigs
      else
         if [[ ${file} =~ ${name_contigs} ]]; then
            # usar ese como ensamble de referencia
            echo "####################################################################################################################"
            echo "  el siguiente genoma se selecciono en base a los valores de N50 de: ${N50_max} y No de contigs de: ${contigs_min}  "
            echo "  se usara como genoma de referencia a = ${file}                                                                    "
            echo "####################################################################################################################"
            echo ""
            # asignar nombre clave de genoma de referencia
            ref_name=$(basename ${file} | cut -d '/' -f '2' | cut -d '.' -f '1' | cut -d '-' -f '1')
            # asignar directorio de salida para guardar resultados de Snippy
            dir="SNIPPY_${ref_name}"
         fi
      fi
   done

   echo "######################################"
   echo -e " Obtener Ensamble de Referencia"
   echo "######################################"

   #################
   # correr Snippy #
   #################

   # for loop para todos los ensambles dentro de la carpeta de ensambles ASSEMBLY
   for ensamble in ASSEMBLY/*.fa; do
      # nombre corto de ensambles
      ensamble_name="$(basename $ensamble .fa)"
      # fila en blanco (espacio)
      echo -e "\n"
      # imprime nombre de ensamble
      echo "#############################################"
      echo "  ${ensamble_name}"
      echo "#############################################"
      echo ""

      # ejecuta Snippy
      snippy --cpus $(nproc) --force --ref ${file} --outdir ${dir}/coreSNP_${ensamble_name} --ctgs ${ensamble} \
      --ram $(grep "MemTotal" /proc/meminfo | awk '{print $2/(1024 * 1024)}' | cut -d "." -f "1") # info de la ram, se divide para Gigas y se eliminan decimales
   done

#####################################################################################
# correr snippy-core, snippy-clean. gubbins, snp-sites, fasttree, raxml y snp-dists #
#####################################################################################

   # moverse al directorio correspondiente
   cd ${dir}

   echo -e "#########################"
   echo " ejecutando snippy core "
   echo -e "#########################\n"

   # ejecuta snippy-core
   snippy-core --ref ../${file} --prefix core $(echo coreSNP_*)

   echo -e "############################################################"
   echo "  limpiando alineamiento para generar SNPs de alta calidad  "
   echo -e "############################################################\n"

   # ejecuta snippy-clean, para limpiar alineamiento
   snippy-clean_full_aln core.full.aln > clean.full.aln
   # correr gubbins con arbol rapido, para limpiar alineamiento
   run_gubbins --verbose --threads $(nproc) --tree_builder fasttree --prefix gubbins clean.full.aln
   # si gubbins falla con error (porque son muy pocas secuencias)
   if [[ $? != 0 ]]; then
      # guardar el resultado de porcentajes de missing data
      run_gubbins --verbose --threads $(nproc) --tree_builder fasttree --prefix gubbins clean.full.aln > tmp_porcentajes.txt
      # filtrar para obtener el valor maximo de porcentajes de missing data encontrados
      maximo=$(cat tmp_porcentajes.txt | grep 'percentage' | awk '{print $NF}' | awk '/[0-9]/' | sort -r | sed -n '1p')
      # volver a correr gubbins con nuevo criterio porcentaje de missing data
      run_gubbins --verbose --threads $(nproc) --filter_percentage $(echo ${maximo} + 1 | bc | cut -d '.' -f '1') --tree_builder fasttree --prefix gubbins clean.full.aln
      # eliminar archivos temporales
      rm tmp*
   fi
   # correr snp-sites, para limpiar alineamiento
   snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln

   echo -e "#############################################"
   echo "  obteniendo reconstrucciones filogeneticas  "
   echo -e "#############################################\n"

   # correr fasttree, para hacer reconstruccion filogenetica por maximum likelihood
   FastTree -gtr -nt clean.core.aln > FastTree_clean.core.tree
   # correr raxml, para hacer segunda reconstruccion filogenetica por maximum likelihood (con 100 bootstraps)
   raxmlHPC -f a -p 1234567890 -s clean.core.aln -x 1234567890 -# 100 -m GTRGAMMA -n clean.core.newick
   # correr snp-dists, para obtener matriz de distancias de SNPs
   snp-dists -j $(nproc) clean.core.aln > Genero_SNP_matrix.tsv
   # limpiar el archivo "Genero_SNP_matrix.tsv" (remover "coreSNP_")
   sed -i 's/coreSNP_//g' Genero_SNP_matrix.tsv

   echo -e "##################################################"
   echo "    Limpiando archivos Newick    "
   echo -e "##################################################"

   # eliminar el string "coreSNP_" del archivo Newick de FastTree
   sed -i 's/coreSNP_//g' FastTree_clean.core.tree
   # eliminar el string "coreSNP_" del archivo Newick de RAxML
   sed -i 's/coreSNP_//g' RAxML_bipartitions.clean.core.newick

   # eliminar el string "-spades-assembly" del archivo Newick de RAxML
   sed -i 's/-spades-assembly//g' RAxML_bipartitions.clean.core.newick
   # eliminar el string "-spades-assembly" del archivo Newick de FastTree
   sed -i 's/-spades-assembly//g' FastTree_clean.core.tree

}

##################################################
# llamar a la funcion que hace todo "run_snippy" #
##################################################
run_snippy

