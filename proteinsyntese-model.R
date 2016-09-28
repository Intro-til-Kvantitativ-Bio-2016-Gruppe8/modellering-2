# Hjælpefunktioner

plotSeveralSeries = function(x, ys) { # Kaldes som følger: plotSeveralSeries(x-værdi-vektor,data.frame(y-værdi-vektor1,y-værdi-vektor2))
  matplot(x,as.matrix(ys), col = 1:length(ys), type="l")
  legend("topright", legend = names(ys), pch = 1, col = 1:length(ys))
}

# Initialisering
t_max = 60*60 # sek i en time
t0 = 1

proteinSyntetiseret = vector(mode="numeric",length=t_max-t0+1)
DproteinSyntetiseret = vector(mode="numeric",length=t_max-t0+1)
translationsProdukter = vector(mode="numeric",length=t_max-t0+1)
maengde_mRNA = vector(mode="numeric",length=t_max-t0+1)
Dmaengde_mRNA = vector(mode="numeric",length=t_max-t0+1)
mRNA_produktion = vector(mode="numeric",length=t_max-t0+1)


# Undermodeller

proteinMasse = function(antalAminosyrer)
{
  antalAminosyrer*118.9
}
avogadro = 6E23
maengde_ribosomer = 3.3E6/avogadro
maengde_tRNA = 12.2*maengde_ribosomer/20 # antager at tRNA er ligeligt fordelt mellem alle typer aminosyrer
ribosomalHastighed = 4 # aminosyrer/sek
andel_aktivt_mRNA = 0.6 # Tennas gæt
mRNA_halfLife = 10*(60*60) # sekunder
mRNA_sekundRente = exp(log(1/2)/mRNA_halfLife) - 1
maengde_aktiv_RNA_polymerase = 60000/avogadro # fra en mus
hastighed_RNA_polymerase_nucleotider = 70 # pr. sekund
gen_længde = 12E3
hastighed_RNA_polymerase_antal_mRNA = hastighed_RNA_polymerase_nucleotider/gen_længde
andel_succesfulde_mRNA = 0.01 # fra BioNumbers


# Model
## Startbetingelser
proteinSyntetiseret[t0] = 0
maengde_mRNA[t0] = 2E5/avogadro
## Simulation
for(t in t0:t_max)
{
  mRNA_produktion[t] = maengde_aktiv_RNA_polymerase*hastighed_RNA_polymerase_antal_mRNA*andel_succesfulde_mRNA
  Dmaengde_mRNA[t] = mRNA_sekundRente*maengde_mRNA[t] + mRNA_produktion[t]
  maengde_mRNA[t+1] = maengde_mRNA[t]+Dmaengde_mRNA[t]
  translationsProdukter[t] = min(andel_aktivt_mRNA*maengde_mRNA[t], maengde_ribosomer, maengde_tRNA)*ribosomalHastighed
  DproteinSyntetiseret[t] = proteinMasse(translationsProdukter[t])
  proteinSyntetiseret[t+1] = proteinSyntetiseret[t] + DproteinSyntetiseret[t]
}
plot(t0:(t_max+1),proteinSyntetiseret,type = "l",
     ylab = "Massen af protein i gram",xlab = "Tiden i sekunder",main = "Massen af protien over tid",
     mgp=c(1.4,.4,0), font.lab=2,
     cex.lab=0.9, cex.axis=0.8,cex.main=1.2)
(10^-9*10^-3)/proteinSyntetiseret[t_max] # Enheden er milliarder af celler der skal til at producere 1 mg protein på en time