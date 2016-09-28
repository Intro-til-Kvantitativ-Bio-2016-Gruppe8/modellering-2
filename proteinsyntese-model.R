# Initialisering
t_max = 60*60 # sek i en time
t0 = 1

proteinSyntetiseret = vector(mode="numeric",length=t_max-t0+1)
DproteinSyntetiseret = vector(mode="numeric",length=t_max-t0+1)
translationsProdukter = vector(mode="numeric",length=t_max-t0+1)
maengde_mRNA = vector(mode="numeric",length=t_max-t0+1)

# Undermodeller

proteinMasse = function(antalAminosyrer)
{
  antalAminosyrer*118.9
}

maengde_ribosomer = 9.98E-18
maengde_tRNA = 12.2*maengde_ribosomer/20 # antager at tRNA er ligeligt fordelt mellem alle typer aminosyrer
ribosomalHastighed = 4 # aminosyrer/sek
andel_aktivt_mRNA = 1

# Model
## Startbetingelser
proteinSyntetiseret[t0] = 0
maengde_mRNA[t0] = 3.321E-19
## Simulation
for(t in t0:t_max)
{
  maengde_mRNA[t+1] = maengde_mRNA[t]
  translationsProdukter[t] = min(andel_aktivt_mRNA*maengde_mRNA[t], maengde_ribosomer, maengde_tRNA)*ribosomalHastighed
  DproteinSyntetiseret[t] = proteinMasse(translationsProdukter[t])
  proteinSyntetiseret[t+1] = proteinSyntetiseret[t] + DproteinSyntetiseret[t]
}
plot(t0:(t_max+1),proteinSyntetiseret,type = "l")
(10^-9*10^-3)/proteinSyntetiseret[t_max] # Enheden er milliarder af celler der skal til at producere 1 mg protein på en time