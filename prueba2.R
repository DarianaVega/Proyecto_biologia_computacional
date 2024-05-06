library(seqinr)


trad =    c(UUU="F", UUC="F", UUA="L", UUG="L",
            UCU="S", UCC="S", UCA="S", UCG="S",
            UAU="Y", UAC="Y", UAA="STOP", UAG="STOP",
            UGU="C", UGC="C", UGA="STOP", UGG="W",
            CUU="L", CUC="L", CUA="L", CUG="L",
            CCU="P", CCC="P", CCA="P", CCG="P",
            CAU="H", CAC="H", CAA="Q", CAG="Q",
            CGU="R", CGC="R", CGA="R", CGG="R",
            AUU="I", AUC="I", AUA="I", AUG="M",
            ACU="T", ACC="T", ACA="T", ACG="T",
            AAU="N", AAC="N", AAA="K", AAG="K",
            AGU="S", AGC="S", AGA="R", AGG="R",
            GUU="V", GUC="V", GUA="V", GUG="V",
            GCU="A", GCC="A", GCA="A", GCG="A",
            GAU="D", GAC="D", GAA="E", GAG="E",
            GGU="G", GGC="G", GGA="G", GGG="G")

ToARN = function(input) {
  newARN = as.vector(input)
  newARN[which(newARN=="t")] = "u"
  newARN = toupper(newARN)
  return (newARN)
}

fRef = read.fasta("wuhan_sequence.txt")
str(fRef)

fMexa = read.fasta("sequence_5.fasta")
str(fMexa)

mutaciones_reg = c()

for (i in seq(1,length(fRef),1)){
  anotaciones = attr(fRef[[i]], "Annot") 
  atributos = unlist(strsplit(anotaciones,"\\[|\\]|:|=|\\.|\\(")); 
  geneName = atributos[which(atributos=="gene")+1] 
  genRef = ToARN( fRef[[i]] )
  for (j in seq(i, length(fMexa),12)) {
    cat(i, "comparado con", j, "\n")
    genMex = ToARN( fMexa[[i]] )  
    cat("Length gen",geneName,":",length(genRef), length(genMex), "\n")
    if (length(genRef) == length(genMex)){
      dif = which(genRef != genMex)
      if (length(dif) > 0){
        for (x in dif){
          if (genRef[x] != "" && genMex[x] !="") {
            muta = paste(genRef[x],"to",genMex[x], sep="") 
            inicioCodon = x - (x-1)%%3 
            numCodon = as.integer((x-1)/3+1) 
            codonOri = paste(genRef[inicioCodon], genRef[inicioCodon+1], genRef[inicioCodon+2],sep="")
            codonMex = paste(genMex[inicioCodon], genMex[inicioCodon+1], genMex[inicioCodon+2],sep="")
            codonChange = paste(codonOri,"to",codonMex, sep="")
            aminoChange = paste(trad[codonOri],numCodon,trad[codonMex], sep="")
            mutacion = paste(muta, codonChange, aminoChange, geneName, sep=" ")
            if (!mutacion %in% mutaciones_reg) {
              cat("\t", mutacion, "\n")
              mutaciones_reg = c(mutaciones_reg, mutacion)
            }
          }
        }
      }
    }else{
      cat("Mutación de deleción o inserción \n")
      genRefLength = length(genRef)
      genMexLength = length(genMex)
      
      if (genRefLength < genMexLength) {
        genRef = c(genRef, rep("", genMexLength - genRefLength))
      } else {
        genMex = c(genMex, rep("", genRefLength - genMexLength))
      }
      dif = which(genRef != genMex)
      if (length(dif) > 0) {
        for (x in dif) {
          if (genRef[x] != "" && genMex[x] !="") {
            muta = paste(genRef[x], "to", genMex[x], sep="") 
            inicioCodon = x - (x - 1) %% 3 
            numCodon = as.integer((x - 1) / 3 + 1) 
            codonOri = paste(genRef[inicioCodon], genRef[inicioCodon + 1], genRef[inicioCodon + 2], sep="")
            codonMex = paste(genMex[inicioCodon], genMex[inicioCodon + 1], genMex[inicioCodon + 2], sep="")
            codonChange = paste(codonOri, "to", codonMex, sep="")
            aminoChange = paste(trad[codonOri], numCodon, trad[codonMex], sep="")
            mutacion = paste(muta, codonChange, aminoChange, geneName, sep=" ")
            if (!mutacion %in% mutaciones_reg) {
              cat("\t", mutacion, "\n")
              mutaciones_reg = c(mutaciones_reg, mutacion)
            }  
          }
        }
      }
    }
  }
}

df = data.frame(
  Mutation = character(),
  Codon = character(),
  Amino = character(),
  Position = integer(),
  Gene = character()
)

obs = list(muta, codonChange, aminoChange, geneName); 
df[10,] = obs
df

library(dplyr)
library(ggplot2)

p = ggplot(df)
p = p + aes(x=Mutation, fill=Mutation, label=after_stat(count))
p = p + ggtitle("Mutaciones de sustitución")
p = p + labs(x="Mutation", y="Count", fill="Count")
p = p + geom_bar(stat = "count")
p = p + geom_text(stat = "count", vjust=1.5)
p = p + facet_grid(~Gene)
p

dfgraph = filter(
  summarise(
    select(
      group_by(df, Amino),
      Mutation:Gene
    ),
    Mutation = first(Mutation),
    Codon = first(Codon),
    Gene = first(Gene),
    Cuenta = n()
  ),
  Cuenta>=1
)

dfgraph

p2 = ggplot(dfgraph)
p2 = p2 + aes(x=Amino, y=Cuenta, fill=Amino, label=Cuenta)
p2 = p2 + ggtitle("Cambio de Aminoácidos")
p2 = p2 + labs(x="Amino", y="Frecuencia", fill="Frecuencia")
p2 = p2 + geom_bar(stat = "identity")
p2 = p2 + geom_text(stat = "identity", vjust=1.5)
p2