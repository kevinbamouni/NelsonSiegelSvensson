# NelsonSiegelSvensson


Calibration du model de Nelson Siegel Svensson à partir des données de la courbe des taux de rendements spot de la BCE
date du 3/3/2017

# Modèle de NELSON-SIEGEL-SVENSSON (NSS)


## Application à la courbes des taux spot de la BCE

Le modèle de NELSON SIEGEL SVENSSON est un modèle d'interpolation de la courbe des taux zéro coupon exprimé en fonction des maturités.


## Présentation du modèle NSS


$$ y(\tau)=\beta_{0}+\beta_{1}\left[\frac{1-exp(\frac{-\tau}{\alpha_{1}})}{\frac{\tau}{\alpha_{1}}}\right]+\beta_{2}\left[\frac{1-exp(\frac{-\tau}{\alpha_{1}})}{\frac{\tau}{\alpha_{1}}}-exp(\frac{-\tau}{\alpha_{1}})\right]+\beta_{4}\left[\frac{1-exp(\frac{-\tau}{\alpha_{2}})}{\frac{\tau}{\alpha_{2}}}-exp(\frac{-\tau}{\alpha_{2}})\right]$$

Ainsi la calibration du modèle se fait par l'estimation des paramètres : $ \beta_{0},\beta_{1},\beta_{2},\beta_{3},\alpha_{1},\alpha_{2}. $ (Avec $ \tau $ la maturité)

Pour estimer les paramètres on résoud un ptobleme d'optimisation qui consiste à minimiser les distances entre les observations théoriques et modélisées: $$min\sum_{k=1}^n (y_{nss}-y_{marché})^{2}$$

