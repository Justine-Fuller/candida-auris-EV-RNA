suppressMessages({library(ggplot2); library(dplyr); library(patchwork)})
setwd("/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV")
E <- "PCG_annotation/Cellular_RNA_analysis/03_Enrichment_results/"
UP <- "#c0392b"; DN <- "#2e6fb7"; NS <- "#c9c8c4"; INK <- "#2b2b2a"; MUT <- "#6b6a66"

## ---- Panel A: volcano of protein-coding genes ----
de <- read.csv("Differential_expression_analysis/03_Cellular_DE_analysis/mRNA_DE_analysis/02_Results/03_mRNA_Cells_DESeq2_Results.csv") %>%
  filter(grepl("^B9J08_", Feature), !is.na(padj))
cap <- 50
de <- de %>% mutate(cls = case_when(padj < 0.05 & log2FoldChange > 0 ~ "Up",
                                    padj < 0.05 & log2FoldChange < 0 ~ "Down", TRUE ~ "NS"),
                    y = pmin(-log10(padj), cap), capped = -log10(padj) > cap)
n_up <- sum(de$cls == "Up"); n_dn <- sum(de$cls == "Down")
cat("tested", nrow(de), "up", n_up, "down", n_dn, "capped", sum(de$capped), "\n")
pA <- ggplot(de, aes(log2FoldChange, y)) +
  geom_point(data = filter(de, cls == "NS"), colour = NS, size = 0.35) +
  geom_point(data = filter(de, cls != "NS"), aes(colour = cls, shape = capped), size = 0.35) +
  scale_colour_manual(values = c(Up = UP, Down = DN), guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
  geom_hline(yintercept = -log10(0.05), linewidth = 0.25, linetype = "dashed", colour = MUT) +
  annotate("text", x = Inf, y = 24, label = paste0("Up\n", format(n_up, big.mark = ",")),
           hjust = 1.05, vjust = 1, size = 2.6, colour = INK, lineheight = 0.9) +
  annotate("text", x = -Inf, y = 24, label = paste0("Down\n", format(n_dn, big.mark = ",")),
           hjust = -0.05, vjust = 1, size = 2.6, colour = INK, lineheight = 0.9) +
  scale_y_continuous(limits = c(0, cap), expand = expansion(mult = c(0, 0.02))) +
  labs(x = expression(log[2]~"fold change"), y = expression(-log[10]~"adjusted"~italic(p))) +
  theme_classic(base_size = 8) +
  theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        axis.text = element_text(colour = MUT), axis.title = element_text(colour = INK))

## ---- Panel B: representative enriched terms ----
rd <- function(f) read.csv(paste0(E, f, ".csv"))
pick <- function(df, ids, dir, src) df %>% filter(ID %in% ids) %>%
  transmute(ID, term = Description, Count, p.adjust, dir, src)
go_up <- rd("ORA_GO_Upregulated"); go_dn <- rd("ORA_GO_Downregulated")
kg_up <- rd("ORA_KEGG_Upregulated"); kg_dn <- rd("ORA_KEGG_Downregulated")
id_of <- function(df, d) { x <- df$ID[df$Description == d]; if (length(x) != 1) stop("term not found: ", d); x }
sel <- bind_rows(
  pick(go_up, c(id_of(go_up, "structural constituent of ribosome"), id_of(go_up, "cytoplasmic translation"),
                id_of(go_up, "amino acid biosynthetic process"), id_of(go_up, "cell wall macromolecule metabolic process")), "Up", "GO"),
  pick(kg_up, c("map03010", "map00190", "map03050", "map01230"), "Up", "KEGG"),
  pick(go_dn, c(id_of(go_dn, "nucleolus"), id_of(go_dn, "RNA processing"), id_of(go_dn, "rRNA processing"),
                id_of(go_dn, "preribosome"), id_of(go_dn, "fungal-type vacuole"),
                id_of(go_dn, "transmembrane transporter activity"), id_of(go_dn, "iron import into cell"),
                id_of(go_dn, "DNA-binding transcription factor activity")), "Down", "GO"),
  pick(kg_dn, c("map03008", "map03040", "map04138"), "Down", "KEGG")) %>%
  mutate(term = sub("Autophagy - yeast", "Autophagy (yeast)", term),
         label = paste0(term, if_else(src == "KEGG", " (KEGG)", "")),
         score = -log10(p.adjust) * if_else(dir == "Up", 1, -1)) %>%
  arrange(dir == "Up", desc(abs(score)))
stopifnot(nrow(sel) == 19)
sel$label <- factor(sel$label, levels = rev(sel$label))
write.csv(sel, "PCG_annotation/Cellular_RNA_analysis/03_Enrichment_results/Fig_caspofungin_transcriptome_terms.csv", row.names = FALSE)
print(sel %>% select(label, Count, p.adjust, dir))
lim <- 35
pB <- ggplot(sel, aes(score, label, fill = dir)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = INK) +
  geom_text(aes(label = Count, hjust = if_else(dir == "Up", -0.25, 1.25)), size = 2.3, colour = MUT) +
  scale_fill_manual(values = c(Up = UP, Down = DN), guide = "none") +
  scale_x_continuous(limits = c(-19, lim), breaks = seq(-15, 30, 5), labels = function(x) abs(x),
                     expand = expansion(mult = 0)) +
  annotate("text", x = -0.5, y = 19.9, label = "Downregulated", hjust = 1, size = 2.6, colour = DN, fontface = "bold") +
  annotate("text", x = 0.5, y = 19.9, label = "Upregulated", hjust = 0, size = 2.6, colour = UP, fontface = "bold") +
  coord_cartesian(clip = "off", ylim = c(0.6, 19.4)) +
  labs(x = expression(-log[10]~"adjusted"~italic(p)), y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.3), axis.ticks.x = element_line(linewidth = 0.3),
        axis.text.y = element_text(colour = INK, size = 7), axis.text.x = element_text(colour = MUT),
        panel.grid.major.x = element_line(colour = "#ecebe8", linewidth = 0.25),
        plot.margin = margin(12, 4, 2, 2))
p <- (pA | pB) + plot_layout(widths = c(1, 1.6)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))
out <- "PCG_annotation/Cellular_RNA_analysis/03_Enrichment_results/Fig_caspofungin_transcriptome"
ggsave(paste0(out, ".png"), p, width = 16, height = 8.5, units = "cm", dpi = 600, bg = "white")
ggsave(paste0(out, ".pdf"), p, width = 16, height = 8.5, units = "cm", bg = "white")
