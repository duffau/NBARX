library(readxl)       # Reading in Excel files
library(quantmod)
library(zoo)          # (Ir)regular spaced tim eseries
library(RColorBrewer) # Color package

# Defining relative paths
# -----------------------
est_output_path <- "../EstOutput/"
# est_output_path <- "est_output/"
data_path <- "../data/"
# data_path <- "data/"
plots_path <- "../plots/"
# plots_path <- "plots/"

# Loading function libraries
# --------------------------
source("func.R")
source("time.utilities.R")

# Defining color vectors
# ----------------------
col_vec <- brewer.pal(9,"Set1")
col_vec2 <- brewer.pal(9,"Pastel1")
col_vec3 <- rgb(c(9,204,164,19,159,0),c(51,0,133,145,159,173),c(83,0,68,35,159,204),maxColorValue=255)
col_vec4 <- rgb(c(0,204,207,19,159,0),c(99,0,184,145,159,173),c(198,0,138,35,159,204),maxColorValue=255)
col_vec5 <- c(col_vec4[1:2],col_vec4[4])

blues <- brewer.pal(9,"Blues")[5:7]
greys <- brewer.pal(9,"Greys")

# Plot output parameters
# ----------------------
w_pdf <- 7
h_pdf <- 5
par_plot <- list(mar=c(2,4,1,1),bg="white")
par_plot_acf <- list(mar=c(4,4,1,1),bg="white")

# Function for page numer in Hamilton (1994) pdf
the_ham <- function(x){
  floor((x+16)/2)
}



