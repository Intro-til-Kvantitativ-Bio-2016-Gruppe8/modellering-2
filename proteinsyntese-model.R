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

maengde_mRNA = 4.3E-10
maengde_ribosomer = 9.98E-18
maengde_tRNA = 6.23E-19
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
proteinSyntetiseret[t_max]