cat("\014")
####################################### CONFIGURACION ####################################### 
# Directorio de los experimentos
dir<-c("groups")
# Significancia estadistica
alpha<-0.05
beta1<-0.10
####################################### CONFIGURACION ####################################### 

## Friedman ranks, aligned ranks en Post Hoc  tests 
## zie : https://github.com/b0rxa/scmamp/blob/master/R/tests.R
## https://cran.r-project.org/web/packages/scmamp/scmamp.pdf 

if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("b0rxa/scmamp")
library("scmamp")

# Nombres de los archivos de cada experimento
files<-list.files(path = dir, pattern = NULL, all.files = FALSE,full.names = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
jf<-length(files)-1

mat<-c()

for (i in 1:length(files)) {
  # Nombre de la prueba
  group<-paste(files[i])
  
  # Obtiene las rutas de cada archivo
  path<-paste(dir,"/",files[i], sep="")
  
  # Lee resultados del experimento ($V1 indica la columna de donde estos se toman)
  value_col <-read.csv(file=path, header=FALSE)$V1
  
  # Crea matriz de datos
  mat <- cbind(mat, value_col)
}
colnames(mat)=files
mat

# Ranking de soluciones por conjunto de datos (por ejecucion) con manejo de empates por promedio
rk<-c()

for( row in 1:nrow(mat)) {
  rk_row<-rank(mat[row,], na.last = TRUE,ties.method = "average")
  rk <- rbind(rk, rk_row)
}

rk <- colSums (rk)/nrow(rk)
rk

# Crea archivo latex
fname <- "post_hoc_friedman.tex"
file.create(fname)
write("\\documentclass{article}\n \\usepackage{float}\n \\usepackage{amsmath} \n \\usepackage[margin=1cm]{geometry}\n \\geometry{a3paper}\n \\begin{document}\n \\title{Post-hoc Friedman tests}\n \\maketitle", file = fname)

# ************************* RANKINGS *****************************************

# Crea tabla para mostrar el ranking
cat(paste("\\begin{table}[tb]\n \\begin{center}
  \\caption{Ranks achieved by the Friedman test in the main case of study. The computed statistics and the related p-value are also shown.}\\label{tb:ranks}
          \\footnotesize
          \\begin{tabular}{ll}
          Algorithms & Ranks \\\\ \\hline \n", sep=""), file = fname, append = TRUE)
for (i in 1:length(files)) {
  cat(paste(files[i]," & ", formatC(rk[i],digits = 4, format = "g"), "\\\\ \n"), file = fname, append = TRUE)
}


ft_result<-friedman.test(mat)
pvalue<-ft_result$p.value
statistic<-ft_result$statistic

#cat(paste("\hline", sep=""), file = fname, append = TRUE)

cat(paste("\\hline\n Statistic & ", formatC(statistic,digits = 4, format = "g"), "\\\\ \n"), file = fname, append = TRUE)
cat(paste("p-value & ", formatC(pvalue,digits = 4, format = "E"), "\\\\ \n"), file = fname, append = TRUE)


pvalue
statistic


# Termina la tabla de ranking
cat("\\end{tabular}\n \\end{center}\n \\end{table}\n", file = fname, append = TRUE)





# ************************* TESTS *****************************************

# Crea una nueva tabla para los resultados
cat(paste("\\begin{table}[tb]\n \\begin{center}
          \\caption{Adjusted p-values for tests for multiple comparisons among all methods.}\\label{tb:",files[i],"}
          \\footnotesize
          \\begin{tabular}{lllll}
          Hypotesis & Unajusted p & Holm & Shaffer & z \\\\ \\hline \n", sep=""), file = fname, append = TRUE)

unadjusted_pval=postHocTest(data=mat, test = "friedman")$raw.pval
holm_corrected_pval=postHocTest(data=mat, test = "friedman", correct = "holm")$corrected.pval
shaffer_corrected_pval=postHocTest(data=mat, test = "friedman", correct = "shaffer")$corrected.pval
# bergmann_corrected_pval=postHocTest(data=mat, test = "friedman", correct = "bergmann")$corrected.pval

k<-ncol(mat)
n<-nrow(mat)
#print("Hypotesis, Unajusted p, Holm, Shaffer, Bergmann, z(Friedman test)")
print("Hypotesis, Unajusted p, Holm, Shaffer, z")
for (i in 1:(length(files)-1)) {
  for (j in (i+1):length(files)) {
    z=(rk[i]-rk[j])/sqrt((k*(k+1))/(6*n))
    #print(paste(files[i]," vs ",files[j], ", ", unadjusted_pval[i,j], ", ", holm_corrected_pval[i,j], ", ", shaffer_corrected_pval[i,j],", ", bergmann_corrected_pval[i,j], ", ", z))
    print(paste(files[i]," vs ",files[j], ", ", unadjusted_pval[i,j], ", ", holm_corrected_pval[i,j], ", ", shaffer_corrected_pval[i,j],", ", z))
    ######## Editado por omar #############################
    if((holm_corrected_pval[i,j]>alpha && holm_corrected_pval[i,j]<beta1) || (shaffer_corrected_pval[i,j]>alpha && shaffer_corrected_pval[i,j]<beta1)){
      if (z<0){
        cat(paste("\\textbf{",files[i],"}  $\\ast \\ast$"," vs  ","\\text{",files[j],"}", " & "), file = fname, append = TRUE)
      } else{
        cat(paste("\\text{",files[i],"}"," vs "," \\textbf{",files[j],"}  $\\ast \\ast$", " & "), file = fname, append = TRUE)
      }
    }
    else if(holm_corrected_pval[i,j] < alpha || shaffer_corrected_pval[i,j] < alpha){
      if (z<0){
        cat(paste("\\textbf{",files[i],"}"," vs  ","\\text{",files[j],"}", " & "), file = fname, append = TRUE)
      } else{
        cat(paste("\\text{",files[i],"}"," vs "," \\textbf{",files[j],"}", " & "), file = fname, append = TRUE)
      }}
    else{
      cat(paste("\\text{",files[i],"}"," vs "," \\text{",files[j],"}", " & "), file = fname, append = TRUE)
    }
    #######################################################
    if(unadjusted_pval[i,j] < alpha) {
      if(unadjusted_pval[i,j]==0||unadjusted_pval[i,j]==1){
        cat(paste("\\textbf{", formatC(unadjusted_pval[i,j],digits = 4, format = "g"), "} & "), file = fname, append = TRUE)
      } else{
        cat(paste("\\textbf{", formatC(unadjusted_pval[i,j],digits = 4, format = "E"), "} & "), file = fname, append = TRUE)
      }
    } else {
      if(unadjusted_pval[i,j]==0||unadjusted_pval[i,j]==1){
        cat(paste("\\text{", formatC(unadjusted_pval[i,j],digits = 4, format = "g"), "} & "), file = fname, append = TRUE)
      }else if ( unadjusted_pval[i,j]>alpha && unadjusted_pval[i,j]<beta1) {
        cat(paste("\\textbf{", formatC(unadjusted_pval[i,j],digits = 4, format = "E"), "} $\\ast \\ast$ & "), file = fname, append = TRUE)
      } 
      else{
        cat(paste("\\text{", formatC(unadjusted_pval[i,j],digits = 4, format = "E"), "} & "), file = fname, append = TRUE)
      }
    }
    #######################################################
    if(holm_corrected_pval[i,j] < alpha) {
      if(holm_corrected_pval[i,j]==0||holm_corrected_pval[i,j]==1){
        cat(paste("\\textbf{", formatC(holm_corrected_pval[i,j],digits = 4, format = "g"), "} & "), file = fname, append = TRUE)
      } else{
        cat(paste("\\textbf{", formatC(holm_corrected_pval[i,j],digits = 4, format = "E"), "} & "), file = fname, append = TRUE)
      }
    } else {
      if(holm_corrected_pval[i,j]==0||holm_corrected_pval[i,j]==1){
        cat(paste("\\text{", formatC(holm_corrected_pval[i,j],digits = 4, format = "g"), "} & "), file = fname, append = TRUE)
      } else if (holm_corrected_pval[i,j]>alpha && holm_corrected_pval[i,j]<beta1) {
        cat(paste("\\textbf{", formatC(holm_corrected_pval[i,j],digits = 4, format = "E"), "} $\\ast \\ast$  & "), file = fname, append = TRUE)
      } 
      else{
        cat(paste("\\text{", formatC(holm_corrected_pval[i,j],digits = 4, format = "E"), "} & "), file = fname, append = TRUE)
      }
    }
    ######################################################
    if(shaffer_corrected_pval[i,j] < alpha) {
      if(shaffer_corrected_pval[i,j]==0||shaffer_corrected_pval[i,j]==1){
        cat(paste("\\textbf{", formatC(shaffer_corrected_pval[i,j],digits = 4, format = "g"), "} & "), file = fname, append = TRUE)
      } else{
        cat(paste("\\textbf{", formatC(shaffer_corrected_pval[i,j],digits = 4, format = "E"), "} & "), file = fname, append = TRUE)
      }
    } else {
      if(shaffer_corrected_pval[i,j]==0||shaffer_corrected_pval[i,j]==1){
        cat(paste("\\text{", formatC(shaffer_corrected_pval[i,j],digits = 4, format = "g"), "} & "), file = fname, append = TRUE)
      } else if (shaffer_corrected_pval[i,j]>alpha && shaffer_corrected_pval[i,j]<beta1) {
        cat(paste("\\textbf{", formatC(shaffer_corrected_pval[i,j],digits = 4, format = "E"), "}  $\\ast \\ast$ & "), file = fname, append = TRUE)
      } 
      else{
        cat(paste("\\text{", formatC(shaffer_corrected_pval[i,j],digits = 4, format = "E"), "} & "), file = fname, append = TRUE)
      }
    }
    #################################################################################
    

    #if(bergmann_corrected_pval[i,j] < alpha) {
     # cat(paste("\\textbf{", bergmann_corrected_pval[i,j], "} & "), file = fname, append = TRUE)
    #} else {
     # cat(paste(bergmann_corrected_pval[i,j], " & "), file = fname, append = TRUE)
    #}
    
    cat(paste( formatC(z,digits = 4, format = "g"), "\\\\ \n"), file = fname, append = TRUE)
  }
}

# Termina la tabla de resultados
cat("\\end{tabular}\n \\end{center}\n \\end{table}\n", file = fname, append = TRUE)

# ********************************************************************************

# Termina archivo latex
cat("\\end{document}", file = fname, append = TRUE)

