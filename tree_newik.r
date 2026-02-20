if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ggtree")
install.packages("ggplot2")
install.packages("ape")

library(ggtree)
library(ggplot2)
library(ape)

# Загрузить дерево после IQ-TREE
tree <- read.tree("supermatrix.fasta.treefile")

#Укоренение дерева 
tree <- root(tree, outgroup = "Outgroup_Sample", resolve.root = TRUE)

# построение
p <- ggtree(tree, layout = "rectangular") + # "circular"-круговое дерево, "fan" - веерное дерево 
    geom_tiplab(size = 4, color = "darkblue", fontface = "italic") + # Названия образцов
    geom_treescale() # Линейка генетической дистанции
#ggtree(tree, branch.length = "none") + geom_tiplab() #для кладограммы

# 4. Добавляем значения Bootstrap (поддержка узлов) выше 70%
p <- p + geom_text2(aes(label = label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),
                    size = 3, color = "red", vjust = -0.5, hjust = 1.1)

# 5. Оформление заголовка и темы
p <- p + ggtitle("Phylogenetic Tree of Hares (BUSCO markers)") +
    theme_tree2() + # Добавляет шкалу внизу
    theme(plot.title = element_text(hjust = 0.5, size = 16))

# 6. Вывод на экран
print(p)

# 7. Сохранение в PDF или PNG высокого качества
ggsave("hare_phylogeny.pdf", p, width = 10, height = 8)
