#source ~/.bashrc

 #my_now.vcf.gz - выравнивание на пищуху, фастк дамп, весь пул
 
 #mY_0502_sort.bam - четыре или сколько там оставшихся до дерева на пищуху


 
busco -i hare_assembly.fasta -o busco_hare_results -l mammalia_odb10 -m genome -c 16
#  -i: ваш входной файл генома (в FASTA).
#-o: название папки, которую создаст BUSCO (не создавать её заранее).
#  -l: база данных (например, mammalia_odb10).
#-m genome: режим поиска в геноме 
#-c 16: количество ядер 

#bed файл
grep -v "^#" ref_busco_krol/run_glires_odb10/full_table.tsv | awk -F '\t' '$2 == "Complete" {print $3"\t"($4 < $5 ? $4 : $5)"\t"($4 < $5 ? $5 : $4)"\t"$1}' > busco_genes.bed

# какие зайцы
bcftools query -l ./send/my_now.vcf.gz > samples.txt

#чета там
./busco.1.py

 conda install -c bioconda mafft
 
#мафт
mkdir -p short_genes
# Копируем только те файлы, которые весят меньше 50 Кб (примерно)
find ortholog_groups/ -name "*.fasta" -size -50k -exec cp {} short_genes/ \;
echo "Отобрано генов: $(ls short_genes | wc -l)"

mkdir -p aligned_short
ls short_genes/*.fasta | xargs -I {} -P 4 sh -c '
    out="aligned_short/$(basename {})"
    mafft --quiet --retree 1 --6merpair {} > "$out"
    echo "Готово: {}"
'

#суперматрица (можно программно)
./busco.3.py

iqtree -s supermatrix.fasta -p partitions.txt -B 1000 -nm 200 -T 10

#ИЛИ
mashtree --sketch-size 10000 --mindepth 0 *.fasta > orthotree.dnd
 #   --mindepth 0: Обязательно добавьте это, если вы подаете на вход уже очищенные FASTA (сборки), а не сырые чтения. Это скажет программе доверять каждому нуклеотиду.
 #   --sketch-size 10000: Сделает «отпечаток» более детальным.
mashtree ortholog_groups/*.fasta > tree.nwk
