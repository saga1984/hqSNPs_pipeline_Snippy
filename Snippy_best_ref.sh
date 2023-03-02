#!/bin/bash

#
# correr snippy en modo de una sola secuencia sobre elementos en una carpeta
#

run_stats(){

# crear directorio de trabajo
mkdir ASSEMBLY/Stats

# volver variable directorio de trabajo
dir="ASSEMBLY/Stats"

# --------------------------------------------------------------------
# primer loop, creacion de archivos de estadisticas con assembly-stats
# --------------------------------------------------------------------
# NOTA: assembly-stats no obtiene profundidad, se debe obtener por separado

for ensamble in ASSEMBLY/*.fa; do
   ename=$(basename ${ensamble} | cut -d "-" -f "1") # asignar nombre de ensamble ( más corto)
   assembly-stats ${ensamble} > ASSEMBLY/Stats/${ename}-stats.txt # crea un archivo de stats por ensamble
done

# ------------------------------------------------------------------------------------
#    segundo loop, creacion del archivo final formato .csv con todas las stadisticas
# ------------------------------------------------------------------------------------

# Generar un archivo y poner los nombres de columnas
echo -e "ID,Contigs,Length,Largest_contig,N50,N90,N_count,Gaps" > ${dir}/estadisticas_ensamble.csv

# ordenar y agrupar estadisticas estadisticas a partir de estadisticas obtenidas por 1er loop
for file in ${dir}/*stats*; do
   # asignar nombre de ensamble (corto)
   fname="$(basename $file | cut -d '-' -f '1')"
   # obtener numero de contigs
   contigs=$(cat ${file} | sed -n '2p' | cut -d ',' -f '2' | cut -d ' ' -f '4')
   # obtener longitud del ensamble
   length=$(cat ${file} | sed -n '2p' | cut -d ',' -f '1' | cut -d ' ' -f '3')
   # obtener tamaño de contig mas grande
   largest=$(cat ${file} | sed -n '2p' | cut -d ',' -f '4' | cut -d ' ' -f '4')
   # obtener N50
   N50=$(cat ${file} | sed -n '3p' | cut -d ',' -f '1' | cut -d ' ' -f '3')
   # obtener N90
   N90=$(cat ${file} | sed -n '7p' | cut -d ',' -f '1' | cut -d ' ' -f '3')
   # obtener Ns
   n_count=$(cat ${file} | sed -n '9p' | cut -d ',' -f '1' | cut -d ' ' -f '3')
   # obtener gaps
   gaps=$(cat ${file} | sed -n '10p' | cut -d ',' -f '1' | cut -d ' ' -f '3')
   # crea filas de archivo con datos obtenidos en las 3 filas anteriores
   echo -e "$fname,$contigs,$length,$largest,$N50,$N90,$n_count,$gaps"
# anexa lo generado por el loop en el archivo creado antes del loop
done >> ${dir}/estadisticas_ensamble.csv

# eliminar '-stats.txt'
sed -i 's/-stats.txt//g' ${dir}/estadisticas_ensamble.csv

# convertir archivo CSV a TSV
sed 's/,/\t/g' ${dir}/estadisticas_ensamble.csv > ${dir}/estadisticas_ensamble.tsv

}

run_snippy() {

   ##########################################
   # elegir el mejor ensamble de referencia #
   ##########################################

   echo -e "\n######################################"
   echo -e " Obtener Ensamble de Referencia"
   echo -e "######################################\n"

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
            echo ""
            echo "####################################################################################################################"
            echo "  el siguiente genoma se selecciono en base a los valores de N50 de: ${N50_max} y No de contigs de: ${contigs_min}  "
            echo "  se usara como genoma de referencia a = ${file}                                                                    "
            echo "####################################################################################################################"
            echo ""
            # asignar nombre clave de genoma de referencia
            ref_name=$(basename ${file} | cut -d '/' -f '2' | cut -d '.' -f '1' | cut -d '-' -f '1')
            # asignar directorio de salida para guardar resultados de Snippy
            dir="SNIPPY_${ref_name}"
            # Asignar como variable al mejor genoma (de referencia)
            Ref=${file}
         fi
      # si los nombres de ensables no son iguales, utiliza solamente el nombre del ensamble obtenido en base a No de contigs
      else
         if [[ ${file} =~ ${name_contigs} ]]; then
            # usar ese como ensamble de referencia
            echo ""
            echo "####################################################################################################################"
            echo "  el siguiente genoma se selecciono en base a los valores de N50 de: ${N50_max} y No de contigs de: ${contigs_min}  "
            echo "  se usara como genoma de referencia a = ${file}                                                                    "
            echo "####################################################################################################################"
            echo ""
            # asignar nombre clave de genoma de referencia
            ref_name=$(basename ${file} | cut -d '/' -f '2' | cut -d '.' -f '1' | cut -d '-' -f '1')
            # asignar directorio de salida para guardar resultados de Snippy
            dir="SNIPPY_${ref_name}"
            # Asignar como variable al mejor genoma (de referencia)
            Ref=${file}
         fi
      fi
   done

   #################
   # correr Snippy #
   #################

   for ensamble in ASSEMBLY/*.fa; do
      # for loop para todos los ensambles dentro de la carpeta de ensambles ASSEMBLY
      ensamble_name="$(basename ${ensamble} -spades-assembly.fa)"

      echo -e "\n####################################################################################################"
      echo -e "  Correr Snippy: para Ensamble = ${ensamble_name}, con Ref = $(basename ${Ref} | cut -d '/' -f '2')  "
      echo -e "####################################################################################################\n"

      # ejecuta Snippy
      snippy --cpus $(nproc) --force --ref ${Ref} --outdir ${dir}/coreSNP_${ensamble_name} --ctgs ${ensamble} \
      --ram $(grep "MemTotal" /proc/meminfo | awk '{print $2/(1024 * 1024)}' | cut -d "." -f "1") # info de la ram, se divide para Gigas y se eliminan decimales
   done
   #####################################################################################
   # correr snippy-core, snippy-clean. gubbins, snp-sites, fasttree, raxml y snp-dists #
   #####################################################################################

   # moverse al directorio correspondiente
   cd ${dir}

   echo -e "\n###################################################################################"
   echo " ejecutando snippy core dentro de: ${dir}, con Ref = $(basename ${Ref} | cut -d '/' -f '2')"
   echo -e "###################################################################################\n"

   # ejecuta snippy-core
   snippy-core --ref ../${Ref} --prefix core $(echo coreSNP_*)

   echo -e "\n############################################################"
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

      echo -e "\n#############################################"
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

      echo -e "\n####################################################"
      echo "   Limpiando archivos Newick y Matriz de distancias   "
      echo -e "####################################################\n"
      # eliminar el string "-spades-assembly" del archivo
      sed -i 's/-spades-assembly//g' Genero_SNP_matrix.tsv

      # eliminar el string "coreSNP_" del archivo Newick de FastTree
      sed -i 's/coreSNP_//g' FastTree_clean.core.tree
      # eliminar el string "coreSNP_" del archivo Newick de RAxML
      sed -i 's/coreSNP_//g' RAxML_bipartitions.clean.core.newick

      # eliminar el string "-spades-assembly" del archivo Newick de RAxML
      sed -i 's/-spades-assembly//g' RAxML_bipartitions.clean.core.newick
      # eliminar el string "-spades-assembly" del archivo Newick de FastTree
      sed -i 's/-spades-assembly//g' FastTree_clean.core.tree
}

###################################################################
# llamar las funciones que hacen todo "run_stats" y  "run_snippy" #
###################################################################
run_stats
run_snippy
