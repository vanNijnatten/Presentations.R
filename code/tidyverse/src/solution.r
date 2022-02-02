library("tidyverse")

data.dir <- "../data"

patient.data <- read.csv(
    file = file.path(data.dir, "patient.data.csv")
  ) %>%
  dplyr::rename(
    patient.id = X
  )
patient.sample.data <- readr::read_csv(
  file = file.path(data.dir, "patient.sample.data.csv")
)
methylation.data <- read.table(
  file = file.path(data.dir, "methylation.data.txt"),
  quote = "",
  sep = "\t"
)
expression.data <- read.table(
  file = file.path(data.dir, "expression.data.txt"),
  sep = "\t"
)
gene.list <- readr::read_csv(
  file = file.path(data.dir, "gene.list.csv")
)

rm("data.dir")


master.table <- patient.data %>%
  dplyr::left_join(
    y = patient.sample.data,
    by = c("patient.id" = "patient.id")
  ) %>%
  dplyr::mutate(
    methylation.id = ifelse(
      paste0("X", patient.id, "_", visit) %in% colnames(methylation.data),
      paste0("X", patient.id, "_", visit),
      NA
    )
  )

methylation.data <- methylation.data %>%
  tibble::rownames_to_column("probe.id") %>%
  tidyr::gather(
    key = "methylation.id",
    value = "methylation.value",
    -probe.id
  )

expression.data <- expression.data %>%
  tibble::rownames_to_column("gene.id") %>%
  tidyr::gather(
    key = "expression.id",
    value = "expression.value",
    -gene.id
  ) %>%
  left_join(
    y = gene.list,
    by = c("gene.id" = "gene.id")
  )

plot.table <- master.table %>%
  dplyr::left_join(
    y = methylation.data,
    by = c("methylation.id" = "methylation.id")
  ) %>%
  dplyr::left_join(
    y = expression.data,
    by = c("expression.id" = "expression.id")
  ) %>%
  dplyr::mutate(
    group = as.factor(group)
  )

gene.to.plot = "gene_1"
probe.to.plot = "probe_1"

plot.table %>%
  dplyr::filter(
    gene.name == gene.to.plot &
    probe.id == probe.to.plot
  ) %>%
  ggplot2::ggplot(
    mapping = aes(
      x = expression.value,
      y = methylation.value,
      color = group
    )
  ) + 
  ggplot2::geom_point()
