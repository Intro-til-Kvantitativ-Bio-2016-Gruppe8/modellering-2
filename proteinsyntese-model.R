# Initialisering
t_max = 60*60 # sek i en time
t0 = 1

proteinSyntetiseret = vector(mode="numeric",length=t_max-t0+1)
translationsProdukter = vector(mode="numeric",length=t_max-t0+1)
proteinSyntese = vector(mode="numeric",length=t_max-t0+1)

# Undermodeller

proteinMasse = function(antalAminosyrer)
{
  antalAminosyrer*118.9
}

maengde_mRNA = 3.321E-19
maengde_ribosomer = 9.98E-18
maengde_tRNA = 12.2*maengde_ribosomer
ribosomalHastighed = 4 # aminosyrer/sek

# Model
## Startbetingelser
proteinSyntetiseret[t0] = 0
## Simulation
for(t in t0:t_max)
{
  translationsProdukter[t] = min(maengde_mRNA, maengde_ribosomer, maengde_tRNA)*ribosomalHastighed
  proteinSyntese[t] = proteinMasse(translationsProdukter[t])
  proteinSyntetiseret[t+1] = proteinSyntetiseret[t] + proteinSyntese[t]
}
plot(t0:(t_max+1),proteinSyntetiseret,type = "l")
(10^-9*10^-3)/proteinSyntetiseret[t_max] # Enheden er milliarder af celler der skal til at producere 1 mg protein på en time