# CHASE

**CHASE** (*CRISPR-Host-phAge Stochastic Evolution*) est un simulateur stochastique écrit en Julia pour étudier la coévolution entre une population bactérienne équipée d'un système CRISPR et une population de phages capables de muter pour y échapper.

C'est un projet personnel que je développe en parallèle de ma formation de niveau M1. L'idée est de partir d'un modèle assez simple pour être compris et testé correctement, puis d'ajouter progressivement du réalisme biologique et des outils d'analyse.

La première version du moteur est fonctionnelle. La partie sur la reconnaissance CRISPR basée sur les séquences est en cours.

## Idée générale

Une bactérie peut conserver dans son locus CRISPR des spacers issus de phages déjà rencontrés. Lorsqu'un phage correspondant réinfecte cette bactérie, il peut être reconnu et neutralisé. À l'inverse, une mutation dans la cible virale peut permettre au phage d'échapper à cette immunité.

On obtient donc une dynamique assez naturelle de course aux armements : de nouveaux clones bactériens apparaissent par acquisition de spacers, tandis que de nouvelles souches phagiques apparaissent par mutation. Ce type de coévolution a déjà été étudié expérimentalement et théoriquement, notamment sur *Streptococcus thermophilus* et le phage 2972 [3,4].

Je ne stocke pas chaque individu séparément. Les bactéries ayant exactement le même contenu CRISPR sont regroupées dans un même clone, et les phages ayant le même génotype dans une même souche.

```julia
Dict{BitSet, Int}  # contenu CRISPR -> effectif du clone
Dict{Int, Int}     # génotype phagique -> effectif de la souche
```

Cela évite de faire grossir inutilement l'état de la simulation lorsque les populations deviennent importantes. En pratique, le coût dépend surtout du nombre de clones et de souches coexistants.

## Événements biologiques

Le modèle contient actuellement six événements :

- division bactérienne ;
- mort bactérienne, avec compétition liée à la capacité de charge ;
- perte d'un spacer ;
- dégradation d'un phage libre ;
- infection reconnue par CRISPR ;
- infection réussie avec lyse de la bactérie et production de nouveaux phages.

Les rencontres entre bactéries et phages suivent une loi d'action de masse : le taux d'interaction entre un clone et une souche est proportionnel au produit de leurs effectifs.

Lors d'une infection réussie, une partie des descendants peut muter et former de nouvelles souches. Une petite probabilité d'acquisition permet également à une bactérie exposée au phage de survivre et d'ajouter un nouveau spacer à son locus CRISPR.

## Moteur de simulation

### Gillespie exact

Le coeur du simulateur utilise le **Stochastic Simulation Algorithm** de Gillespie [1].

À chaque itération, le moteur calcule les taux des événements possibles, tire le temps jusqu'au prochain événement dans une loi exponentielle, choisit un événement proportionnellement à son taux puis met à jour l'état.

L'intérêt principal est de conserver une vraie dynamique stochastique, surtout lorsque les populations deviennent faibles. C'est particulièrement important près des extinctions, où une approximation déterministe ferait facilement disparaître des effets de hasard qui comptent réellement.

### `RatesCache`

Le principal problème de performance vient du nombre de couples clone/souche possibles.

Un recalcul complet de tous les taux après chaque événement devient rapidement coûteux lorsque la diversité augmente. `RatesCache` conserve donc les valeurs déjà calculées et ne met à jour que celles qui sont réellement affectées par le dernier événement.

Le taux total est corrigé par différence au lieu d'être recalculé entièrement à chaque étape. Une resynchronisation complète est faite périodiquement pour limiter l'accumulation d'erreurs d'arrondi.

Le tirage dans `sample_event` a également été écrit pour éviter les allocations inutiles dans la boucle principale.

### Tau-leaping

Lorsque les populations deviennent très grandes, Gillespie exact devient beaucoup moins intéressant : le moteur passe énormément de temps à simuler des événements individuels qui modifient très peu l'état global.

Dans ce cas, CHASE peut temporairement passer en **tau-leaping** [2]. Sur un petit intervalle de temps, plusieurs occurrences d'une même réaction sont tirées en une seule fois avec une loi de Poisson.

Le pas est volontairement limité pour éviter des mises à jour trop brutales. J'ai aussi ajouté des garde-fous pour empêcher qu'un ensemble de réactions tirées pendant le même pas consomme plus d'individus que ceux réellement disponibles.

Lorsque les populations redescendent, le simulateur repasse en Gillespie exact.

## État actuel

La couche de base est fonctionnelle et les tests actuels passent.

Une simulation complète prend typiquement **40 à 60 secondes** sur la machine utilisée pendant le développement, selon les paramètres et surtout selon la diversité produite.

Les simulations montrent déjà des successions d'acquisition de spacers et de mutations d'échappement qui ressemblent qualitativement à une dynamique de type Red Queen. Je considère cependant cela uniquement comme une validation du comportement du modèle. La comparaison à des données biologiques réelles viendra plus tard.

## Structure du projet

```text
src/
└── core/
    ├── state.jl
    ├── events.jl
    ├── rates.jl
    └── gillespie.jl

scripts/
└── run_simulation.jl

test/
└── runtests.jl
```

- `state.jl` contient la représentation des clones, des souches et de la reconnaissance.
- `events.jl` applique les événements biologiques.
- `rates.jl` contient le calcul des taux et le cache incrémental.
- `gillespie.jl` contient la boucle principale et le passage éventuel en tau-leaping.
- `run_simulation.jl` sert à lancer une simulation et générer les figures.
- `runtests.jl` regroupe les tests du moteur.

## Lancer une simulation

Installer les dépendances :

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Puis lancer :

```bash
julia --project=. scripts/run_simulation.jl
```

Le script génère actuellement `recap.png`, avec les trajectoires bactériennes et phagiques ainsi que l'évolution du nombre de clones et de souches.

Pour lancer les tests :

```bash
julia --project=. test/runtests.jl
```

## Roadmap

Je préfère développer le projet par couches plutôt que d'ajouter tout de suite beaucoup de mécanismes biologiques. L'idée est que chaque étape reste testable et compréhensible avant de passer à la suivante.

### Couche 0 - Simulateur Gillespie ✅

Cette partie correspond au moteur actuel.

Elle comprend Gillespie exact, le tau-leaping hybride, les six événements biologiques de base, `RatesCache`, les `BitSet` pour représenter les loci CRISPR et le tirage sans allocation dans la boucle principale.

Le modèle reste volontairement abstrait. La reconnaissance est encore binaire et les génotypes phagiques ne correspondent pas encore à de vraies séquences.

Ce choix est assez proche des premiers modèles de dynamique de populations CRISPR-phage, qui cherchent d'abord à comprendre les grandes dynamiques de coexistence, de résistance et d'échappement avant de représenter toute la biologie moléculaire [3].

### Couche 1 - Reconnaissance CRISPR basée sur les séquences 🔧

C'est la partie en cours.

Le premier changement est de supprimer progressivement la reconnaissance purement binaire. Les phages auront de petites séquences représentant leurs protospacers et les bactéries auront les séquences correspondant à leurs spacers.

Une mutation phagique modifiera donc une séquence existante au lieu de créer simplement un nouvel identifiant.

Dans une première version, la probabilité de reconnaissance dépendra du nombre de mismatches et de leur position. Je ne veux pas utiliser une simple distance de Hamming comme vérité biologique : les données montrent clairement que tous les mismatches n'ont pas le même effet, notamment selon leur position par rapport au PAM [5,6].

Le but est donc surtout d'avoir un modèle intermédiaire : plus réaliste que le tout-ou-rien actuel, mais encore assez simple pour rester rapide et interprétable.

Le PAM pourra ensuite être ajouté explicitement.

### Couche 2 - Inférence statistique 🔲

Une fois le simulateur suffisamment stable, je veux essayer de l'utiliser dans l'autre sens.

Au lieu de fixer les paramètres puis de regarder ce qu'il se passe, l'idée sera de partir d'une trajectoire observée et d'estimer quels paramètres ont pu la produire.

Les paramètres concernés pourront inclure par exemple les taux de naissance et de mort, le taux d'adsorption, le taux de mutation, la probabilité d'acquisition d'un spacer, la perte de spacers ou encore la taille du burst phagique.

Comme la vraisemblance exacte du modèle devient vite difficile à manipuler, une première approche sera l'**Approximate Bayesian Computation** [7,8].

Le principe sera de lancer beaucoup de simulations avec des paramètres différents puis de comparer des statistiques simples : densité moyenne, temps avant extinction, amplitude des oscillations, diversité des clones, diversité des souches ou vitesse d'apparition de nouveaux génotypes.

Ce sera probablement la première partie du projet où le coût de calcul deviendra vraiment limitant.

### Couche 3 - Surrogate model et simulation-based inference 🔲

Cette couche dépend directement de la précédente.

Si une simulation prend plusieurs dizaines de secondes, faire de l'inférence avec des dizaines de milliers de trajectoires devient vite trop coûteux. Je veux donc tester deux pistes.

La première est un **surrogate model** capable d'imiter certaines sorties de CHASE beaucoup plus rapidement. Une Neural ODE est une possibilité intéressante pour apprendre des trajectoires agrégées [9], même si elle ne remplacera probablement pas toute la stochasticité du moteur.

La seconde piste est la **simulation-based inference**. L'idée est d'entraîner un modèle à apprendre directement la relation entre les paramètres du simulateur et les observations, au lieu de rejeter une grande partie des simulations comme en ABC.

Les normalizing flows sont particulièrement intéressants pour approximer des distributions de paramètres dans ce type de problème [10].

Je considère cette partie comme exploratoire : elle n'aura de sens qu'une fois le simulateur biologique suffisamment stable.

### Couche 4 - Structure spatiale 🔲

Le modèle actuel suppose un milieu parfaitement mélangé.

À terme, je veux pouvoir séparer la population en plusieurs patches connectés par migration. Chaque patch aura ses propres bactéries et phages, avec des infections locales, tandis qu'une fraction des individus pourra passer d'un patch à l'autre.

Cela permettra de regarder des phénomènes absents du modèle actuel : extinctions locales, recolonisations, refuges bactériens, vagues d'invasion et maintien de diversité entre zones.

Cette extension est directement motivée par la littérature sur la coexistence phage-bactérie en environnement spatial. Des modèles stochastiques ont montré que des refuges spatiaux peuvent stabiliser la coexistence [11], et des expériences plus récentes montrent que migration et structure spatiale peuvent maintenir plusieurs étapes successives de diversification [12].

### Couche 5 - Séquences et données réelles 🔲

Les séquences utilisées dans la Couche 1 resteront d'abord simplifiées.

Une étape suivante sera d'utiliser de vraies séquences nucléotidiques. **CRISPRCasdb** pourra servir de source pour les arrays CRISPR et les gènes `cas` [13], tandis que **IMG/VR** fournit une grande collection de génomes et fragments de génomes viraux [14].

Le matching pourra alors prendre en compte de vrais protospacers, les mismatches et le PAM.

Un système particulièrement intéressant pour tester cette partie est celui étudié par Common et al. sur *Streptococcus thermophilus* DGCC7710 et le phage 2972 [4]. Les auteurs ont suivi les populations pendant plusieurs semaines et observé à la fois l'acquisition de nouveaux spacers et des mutations dans les protospacers et les PAM des phages.

Ce serait un bon premier cas pour comparer CHASE à autre chose qu'à ses propres trajectoires.

L'objectif initial ne sera pas de reproduire chaque courbe expérimentale point par point. Je veux d'abord voir si le modèle retrouve les mêmes tendances générales : acquisition de résistance, échappement phagique, maintien ou perte de diversité et extinction éventuelle d'une des deux populations.

### Couche 6 - Anti-CRISPR 🔲

La dernière extension prévue pour l'instant concerne les protéines anti-CRISPR.

Certains phages peuvent produire des protéines capables d'inhiber directement la machinerie CRISPR-Cas [15]. Ce mécanisme est différent d'une simple mutation d'échappement, puisqu'un phage peut rester reconnaissable au niveau de sa séquence tout en perturbant la réponse immunitaire.

Les anti-CRISPR sont aussi intéressants d'un point de vue de dynamique de population. Landsberger et al. ont montré expérimentalement que plusieurs phages Acr+ peuvent coopérer : une première infection peut affaiblir la défense CRISPR d'une bactérie et faciliter une infection suivante [16].

Dans CHASE, cela permettrait donc d'ajouter une vraie stratégie phagique supplémentaire, à côté des mutations de protospacers.

À terme, le modèle pourrait comparer des souches qui échappent principalement par mutation, des souches portant un système anti-CRISPR, et éventuellement des souches combinant les deux.

## Questions que j'aimerais explorer

À terme, le projet devrait surtout servir à poser des questions quantitatives :

- dans quelles conditions bactéries et phages coexistent-ils durablement ?
- quand obtient-on plutôt une extinction rapide ?
- quel est l'effet d'un coût de fitness associé au système CRISPR ?
- quelle quantité de diversité est maintenue par la structure spatiale ?
- quels paramètres peuvent réellement être retrouvés à partir d'une trajectoire observée ?
- est-ce qu'un surrogate model peut accélérer l'inférence sans perdre les extinctions rares ?
- dans quelles conditions un anti-CRISPR devient-il avantageux par rapport à une stratégie basée uniquement sur les mutations d'échappement ?

## Références

1. Gillespie, D. T. (1977). **Exact stochastic simulation of coupled chemical reactions.** *The Journal of Physical Chemistry*, 81(25), 2340-2361. https://doi.org/10.1021/j100540a008

2. Gillespie, D. T. (2001). **Approximate accelerated stochastic simulation of chemically reacting systems.** *The Journal of Chemical Physics*, 115, 1716-1733. https://doi.org/10.1063/1.1378322

3. Levin, B. R., Moineau, S., Bushman, M. & Barrangou, R. (2013). **The Population and Evolutionary Dynamics of Phage and Bacteria with CRISPR-Mediated Immunity.** *PLOS Genetics*, 9(3), e1003312. https://doi.org/10.1371/journal.pgen.1003312

4. Common, J., Morley, D., Westra, E. R. & van Houte, S. (2019). **CRISPR-Cas immunity leads to a coevolutionary arms race between Streptococcus thermophilus and lytic phage.** *Philosophical Transactions of the Royal Society B*, 374, 20180098. https://doi.org/10.1098/rstb.2018.0098

5. Farasat, I. & Salis, H. M. (2016). **A Biophysical Model of CRISPR/Cas9 Activity for Rational Design of Genome Editing and Gene Regulation.** *PLOS Computational Biology*, 12(1), e1004724. https://doi.org/10.1371/journal.pcbi.1004724

6. Singh, D. et al. (2016). **Real-time observation of DNA recognition and rejection by the RNA-guided endonuclease Cas9.** *Nature Communications*, 7, 12778. https://doi.org/10.1038/ncomms12778

7. Beaumont, M. A., Zhang, W. & Balding, D. J. (2002). **Approximate Bayesian Computation in Population Genetics.** *Genetics*, 162(4), 2025-2035. https://doi.org/10.1093/genetics/162.4.2025

8. Warne, D. J., Baker, R. E. & Simpson, M. J. (2019). **Simulation and inference algorithms for stochastic biochemical reaction networks: from basic concepts to state-of-the-art.** *Journal of the Royal Society Interface*, 16, 20180943. https://doi.org/10.1098/rsif.2018.0943

9. Chen, R. T. Q., Rubanova, Y., Bettencourt, J. & Duvenaud, D. (2018). **Neural Ordinary Differential Equations.** *Advances in Neural Information Processing Systems 31*.

10. Tejero-Cantero, A. et al. (2020). **sbi: A toolkit for simulation-based inference.** *Journal of Open Source Software*, 5(52), 2505. https://doi.org/10.21105/joss.02505

11. Heilmann, S., Sneppen, K. & Krishna, S. (2012). **Coexistence of phage and bacteria on the boundary of self-organized refuges.** *Proceedings of the National Academy of Sciences*, 109(31), 12828-12833. https://doi.org/10.1073/pnas.1200771109

12. Shaer Tamar, E. & Kishony, R. (2022). **Multistep diversification in spatiotemporal bacterial-phage coevolution.** *Nature Communications*, 13, 7971. https://doi.org/10.1038/s41467-022-35351-w

13. Pourcel, C. et al. (2020). **CRISPRCasdb: a successor of CRISPRdb containing CRISPR arrays and cas genes from complete genome sequences.** *Nucleic Acids Research*, 48(D1), D535-D544. https://doi.org/10.1093/nar/gkz915

14. Camargo, A. P. et al. (2023). **IMG/VR v4: an expanded database of uncultivated virus genomes within a framework of extensive functional, taxonomic, and ecological metadata.** *Nucleic Acids Research*, 51(D1), D733-D743. https://doi.org/10.1093/nar/gkac1037

15. Bondy-Denomy, J., Pawluk, A., Maxwell, K. L. & Davidson, A. R. (2013). **Bacteriophage genes that inactivate the CRISPR/Cas bacterial immune system.** *Nature*, 493, 429-432. https://doi.org/10.1038/nature11723

16. Landsberger, M. et al. (2018). **Anti-CRISPR Phages Cooperate to Overcome CRISPR-Cas Immunity.** *Cell*, 174(4), 908-916.e12. https://doi.org/10.1016/j.cell.2018.05.058

## Dépendances

Le projet est écrit en Julia. La principale dépendance utilisée actuellement est [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl), pour les lois exponentielle, binomiale et de Poisson.

Le projet est encore en développement. La roadmap décrit la direction actuelle du travail et pourra évoluer à mesure que le modèle est testé.
