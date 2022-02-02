library("tidyverse")

set.seed(1)

# Data
patient.data <- tibble::tibble(
  patient.id = 1:10,
  age = c(40, 45, 50, 55, 60, 65, 50, 55, 60, 65),
  group = c(
    rep("treatment", 5),
    rep("plaebo", 5)
  )
)


patient.sample.data <- tibble::tibble(
  patient.id = rep(1:10, 3), 
  visit = c(
    rep(1, 10),
    rep(2, 10),
    rep(3, 10)
  ), 
  treatment = c(
    rep("none", 10),
    rep("medicine", 15),
    rep("placebo", 5)
  )
)

patient.sample.data <- patient.sample.data %>%
  dplyr::left_join(
    y = patient.sample.data %>%
      dplyr::select(
        patient.id,
        visit
      ) %>%
      dplyr::mutate(
        expression.id = paste0(
          "expr_",
          seq(from = 55, to = 144, by = 3)
        )
      ) %>%
      dplyr::sample_n(
        size = 25,
        replace = FALSE
      ),
    by = c(
      "patient.id" = "patient.id",
      "visit" = "visit"
    )
  )
  


rownames <- paste0("probe_", 1:1000)
colnames <- patient.sample.data %>%
  dplyr::sample_n(
    size = 25,
    replace = FALSE
  ) %>%
  dplyr::mutate(
  	sample.id =
  	  paste0("X", patient.id, "_", visit)
  ) %>%
  dplyr::pull(sample.id)

methylation.data <- matrix(
  rnorm(length(colnames) * length(rownames)),
  nrow = length(rownames),
  byrow = TRUE,
  dimnames = list(
    rownames,
    colnames
  )
)
rm(list = c("rownames", "colnames"))


rownames <- paste0(1:500, "_st")
colnames <- patient.sample.data %>% dplyr::filter(!is.na(expression.id)) %>% pull(expression.id)

expression.data <- matrix(
  rnorm(length(colnames) * length(rownames)),
  nrow = length(rownames),
  byrow = TRUE,
  dimnames = list(
    rownames,
    colnames
  )
)
rm(list = c("rownames", "colnames"))


gene.list <- data.frame(
  gene.id = row.names(expression.data),
  gene.name = paste0("gene_", 1:nrow(expression.data))
)


# Write data for later use
data.dir <- "../data"
dir.create(data.dir)

write.csv(
  x = patient.data %>% tibble::column_to_rownames("patient.id"),
  file = file.path(data.dir, "patient.data.csv")
)
readr::write_csv(
  x = patient.sample.data,
  path = file.path(data.dir, "patient.sample.data.csv")
)
write.table(
  x = methylation.data %>% as.data.frame(),
  file = file.path(data.dir, "methylation.data.txt"),
  quote = FALSE,
  sep = "\t"
)
write.table(
  x = expression.data %>% as.data.frame(),
  file = file.path(data.dir, "expression.data.txt"),
  sep = "\t"
)
readr::write_csv(
  x = gene.list,
  path = file.path(data.dir, "gene.list.csv")
)

rm("data.dir")