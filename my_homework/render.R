cmd.args <- commandArgs(trailingOnly = TRUE)

for(f in cmd.args){
    rmarkdown::render(f)
    f.rmd <- gsub("R$", "Rmd", f)
    file.remove(f.rmd)
}
