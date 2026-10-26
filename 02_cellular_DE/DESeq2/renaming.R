library(tidyverse)

# Only Confirmed/High confidence mappings from your JBrowse checks.
# Anything Medium confidence or lower is intentionally left out — those IDs
# will pass through unchanged rather than risk a wrong rename.
id_mapping <- tibble::tribble(
  ~old_name,       ~correct_name,
  "B9J08_00417",   "B9J08_000431",
  "B9J08_00418",   "B9J08_000431",
  "B9J08_00420",   "B9J08_000432",
  "B9J08_00423",   "B9J08_000435",
  "B9J08_00424",   "B9J08_000435",
  "B9J08_00425",   "B9J08_000435",
  "B9J08_00428",   "B9J08_000436",
  "B9J08_00435",   "B9J08_000441",
  "B9J08_00436",   "B9J08_000441",
  "B9J08_00437",   "B9J08_000442",
  "B9J08_00440",   "B9J08_000444",
  "B9J08_00443",   "B9J08_000446",
  "B9J08_00444",   "B9J08_000446",
  "B9J08_03310",   "B9J08_000315",
  "B9J08_03311",   "B9J08_000315"
)

setwd("/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline")
# Read your existing down-regulated gene ID list
gene_ids <- read_lines("Gene_IDs_Downregulated.txt")

renamed_ids <- tibble(old_name = gene_ids) %>%
  left_join(id_mapping, by = "old_name") %>%
  mutate(final_name = coalesce(correct_name, old_name)) %>%
  pull(final_name)

write_lines(renamed_ids, "Gene_IDs_Downregulated_renamed_Annotated.txt")

cat("Renamed:", sum(gene_ids %in% id_mapping$old_name), "of", nrow(id_mapping), "possible IDs\n")
cat("Unchanged (not in mapping table, includes Medium/unconfirmed/real gaps):",
    sum(!gene_ids %in% id_mapping$old_name), "\n")
# checking for duplicates
renamed_ids[renamed_ids == "B9J08_000431"]
sum(duplicated(renamed_ids))  
#collapsing duplicates
renamed_ids_deduped <- unique(renamed_ids)

cat("Before dedup:", length(renamed_ids), "\n")
cat("After dedup:", length(renamed_ids_deduped), "\n")

write_lines(renamed_ids_deduped, "Gene_IDs_Downregulated_renamed_Annotated.txt")

#save
write_lines(
  renamed_ids_deduped,
  "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/Gene_IDs_Downregulated_renamed_Annotated.txt"
)

# confirming where it went
file.exists("/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/Gene_IDs_Downregulated_renamed_Annotated.txt")



# Same Confirmed/High confidence mapping table as the downregulated version now for up
id_mapping <- tibble::tribble(
  ~old_name,       ~correct_name,
  "B9J08_00417",   "B9J08_000431",
  "B9J08_00418",   "B9J08_000431",
  "B9J08_00420",   "B9J08_000432",
  "B9J08_00423",   "B9J08_000435",
  "B9J08_00424",   "B9J08_000435",
  "B9J08_00425",   "B9J08_000435",
  "B9J08_00428",   "B9J08_000436",
  "B9J08_00435",   "B9J08_000441",
  "B9J08_00436",   "B9J08_000441",
  "B9J08_00437",   "B9J08_000442",
  "B9J08_00440",   "B9J08_000444",
  "B9J08_00443",   "B9J08_000446",
  "B9J08_00444",   "B9J08_000446",
  "B9J08_03310",   "B9J08_000315",
  "B9J08_03311",   "B9J08_000315"
)

# Read your UPregulated gene ID list this time
up_gene_ids <- read_lines("Gene_IDs_Upregulated.txt")

renamed_up_ids <- tibble(old_name = up_gene_ids) %>%
  left_join(id_mapping, by = "old_name") %>%
  mutate(final_name = coalesce(correct_name, old_name)) %>%
  pull(final_name)

cat("Before dedup:", length(renamed_up_ids), "\n")
renamed_up_ids_deduped <- unique(renamed_up_ids)
cat("After dedup:", length(renamed_up_ids_deduped), "\n")

write_lines(
  renamed_up_ids_deduped,
  "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/Gene_IDs_Upregulated_renamed_Annotated.txt"
)

file.exists("/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/Gene_IDs_Upregulated_renamed_Annotated.txt")
length(read_lines("/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/Gene_IDs_Upregulated_renamed_Annotated.txt"))
