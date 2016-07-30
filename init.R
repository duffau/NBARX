library(readxl)       # Reading in Excel files
# library(quantmod)
library(zoo)          # (Ir)regular spaced tim eseries
library(RColorBrewer) # Color package
setwd("~/Dropbox/Speciale/r_scripts")


# Defining relative paths
# -----------------------
# est_output_path <- "../EstOutput/"
est_output_path <- "est_output/"
# data_path <- "../data/"
data_path <- "data/"
# plots_path <- "../plots/"
plots_path <- "plots/"

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
par_plot <- list(mar=c(2,4,1,1),bg="white",cex=1.1)
par_plot_acf <- list(mar=c(4,4,1,1),bg="white",cex=1.1)

# Thesis specific utility functions
# ---------------------------------

# Function for page numer in Hamilton (1994) pdf
the_ham <- function(x){
  floor((x+16)/2)
}

# Format columnnames in LaTeX tables
format_colnames <- function(col_names){
  k <- length(col_names)
  colnm <-  rep(NA, k)
  for(i in 1:k){
    str <- col_names[i]
    id1 <- 0
    for(check_str in c("omega","alpha","beta","nu")) id1 <- id1 + sum(grepl(check_str,str))
    if(id1>0){
      id2 <- 0
      for(check_str in c("alpha","beta")) id2 <- id2 + sum(grepl(check_str,str))
      if(id2>0){
        str <- gsub("_", "_{", str, fixed = TRUE)
        colnm[i] <- paste0("$\\",str,"}$")
      } else {
        colnm[i] <- paste0("$\\",str,"$")
      }
    } else {
      str <- gsub("_", "\\_", str, fixed = TRUE)
      colnm[i] <- paste0("$\\mathtt{",str,"}$")
    }
  }
  return(colnm)
}


