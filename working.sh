
 my_now.vcf.gz - выравнивание на пищуху, фастк дамп, весь пул
 
 mY_0502_sort.bam - четыре или сколько там оставшихся до дерева на пищуху

 conda install -c bioconda mafft
 
 
busco -i hare_assembly.fasta -o busco_hare_results -l mammalia_odb10 -m genome -c 16
\
Разбор параметров:

    -i: ваш входной файл генома (в FASTA).

    -o: название папки, которую создаст BUSCO (не создавайте её заранее!).

    -l: база данных (например, mammalia_odb10).

    -m genome: режим поиска в геноме (находит гены «с нуля»).№№№ есть автоматичн

    -c 16: количество ядер (на сервере можно ставить 16, 32 или больше, чтобы было быстрее).

grep -v "^#" ref_busco_krol/run_glires_odb10/full_table.tsv | awk -F '\t' '$2 == "Complete" {print $3"\t"($4 < $5 ? $4 : $5)"\t"($4 < $5 ? $5 : $4)"\t"$1}' > busco_genes.bed

reference.fasta (индексирован: samtools faidx reference.fasta).

    your_data.vcf.gz (индексирован: bcftools index your_data.vcf.gz).

    samples.txt (список имен зайцев из VCF).

    busco_genes.bed (получен из full_table.tsv).

    Важный момент по формату: Если в вашем BED-файле колонок больше (например, 6), а в read вы написали только 4, то всё, что идет после 4-й колонки, Bash «свалит» в последнюю переменную (gene). Для нашей задачи это не страшно, но если колонок меньше, скрипт может работать неправильно.

Хотите проверить, сколько реально колонок в вашем BED-файле, чтобы мы убедились, что скрипт их правильно прочитает? Можете прислать первые две-три строки из него (команда head -n 3 busco_genes.bed).
bcftools query -l ./send/my_now.vcf.gz > samples.txt

#!/bin/bash

# --- НАСТРОЙКИ ---
REF="reference.fasta"
VCF="your_data.vcf.gz"
BED="busco_genes.bed"
SAMPLES="samples.txt"
THREADS=8

# --- ЭТАП 1: Извлечение ортологов из VCF ---
echo "--- Шаг 1: Экстракция генов ---"
mkdir -p ortholog_groups

while read sample; do
    echo "Обработка особи: $sample"
    while read chrom start end gene; do
        # Пишем заголовок (имя особи)
        echo ">${sample}" >> "ortholog_groups/${gene}.fasta"
        # Вырезаем ген, накладываем мутации из VCF, берем только строку букв
        samtools faidx "$REF" "${chrom}:${start}-${end}" | \
        bcftools consensus -s "$sample" -e 'N' "$VCF" | \
        grep -v "^>" >> "ortholog_groups/${gene}.fasta"
    done < "$BED"
done < "$SAMPLES"

# --- ЭТАП 2: Множественное выравнивание (MAFFT) ---
echo "--- Шаг 2: Выравнивание (MSA) ---"
mkdir -p aligned_genes
for f in ortholog_groups/*.fasta; do
    gene_name=$(basename "$f")
    # Проверяем, чтобы файл не был пустым
    if [ -s "$f" ]; then
        mafft --auto --thread "$THREADS" "$f" > "aligned_genes/${gene_name}"
    fi
done

# --- ЭТАП 3: Склейка в суперматрицу (AMAS) ---
echo "--- Шаг 3: Склейка генов ---"
# Установка AMAS: pip install AMAS
python3 -m AMAS concat -i aligned_genes/*.fasta -f fasta -d dna -u fasta -t concatenated.fasta -p partitions.txt

# --- ЭТАП 4: Построение дерева (RAxML-NG) ---
echo "--- Шаг 4: Построение дерева ---"
# Мы используем модель GTR+G (стандарт для ДНК) и 100 реплик бутстрепа
raxml-ng --all --msa concatenated.fasta --model GTR+G --prefix my_hare_tree --threads "$THREADS" --bs-trees 100
