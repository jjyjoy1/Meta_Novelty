library(limma)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Load and prepare data
data <- read.csv(input_file, row.names = 1)
group <- factor(c("Control", "Treatment", "Treatment"))  # Customize as needed

# Create design matrix
design <- model.matrix(~ group)
fit <- lmFit(as.matrix(data), design)
fit <- eBayes(fit)

# Extract results
results <- topTable(fit, coef=2, number=Inf, adjust="fdr")

# Save
write.csv(results, file=output_file)


