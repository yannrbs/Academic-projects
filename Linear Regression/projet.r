setwd("/Users/yann/Documents/Etude sup/INSA MA 3A/Cours/S6/Regression linéaire/Projet")

# Librairies nécessaires pour le test de Durbin-Watson et VIF
library(lmtest)

# Chargement des données
data_base = read.csv2("3MA_Projets_RegressionLinaire.csv")

# Sélection des 20 premières années (données historiques connues)
data = data_base[1:20,]

# Passage au logarithme pour stabiliser les grandes valeurs
data_log = data
data_log[,-1] = log(data[,-1])
summary(data_log)

# Statistiques descriptives et corrélations
cor(data_log)
pairs(data_log)

# Première recherche : modèle complet sur toutes les années
# Élimination pas à pas des variables qui ont leur p-value > 0.05
log_reg1 = lm(TotalPassengers_US~ . - Annee, data_log)
summary(log_reg1)

# On enlève CPI (plus grande p-value)
log_reg2 = lm(TotalPassengers_US~ . - Annee - CPI, data_log)
summary(log_reg2)

# On enlève Conso_gouv
log_reg3 = lm(TotalPassengers_US~ . - Annee - CPI - Conso_gouv, data_log)
summary(log_reg3)

# On enlève GDP
log_reg4 = lm(TotalPassengers_US~ . - Annee - CPI - Conso_gouv - GDP, data_log)
summary(log_reg4)

# On enlève Exports 
log_reg5 = lm(TotalPassengers_US~ . - Annee - CPI - Conso_gouv - GDP - Exports, data_log)
summary(log_reg5)

# On enlève Conso_Privee pour obtenir le modèle réduit
log_reg6 = lm(TotalPassengers_US~ . - Annee - CPI - Conso_gouv - GDP - Exports - Conso_Privee, data_log)
summary(log_reg6)

# Test de Fisher pour valider la réduction du modèle
anova(log_reg1, log_reg6)
# La p-value > 0.05 valide notre modèle réduit

# Analyse des résidus
residu_log = log_reg6$residuals
sum(residu_log)
# La somme des résidus est bien égale à 0

# Test d'autocorrélation de Durbin-Watson
dwtest(log_reg6)
# On a bien p-value > 0.05 donc on ne rejette pas H0
# on a quand meme p-value = 0.08, proche de alpha 

# Affichage des graphiques de diagnostic de base
x11()
par(mfrow = c(2,2))
plot(log_reg6)
par(mfrow = c(1,1))

# Diagnostic des points leviers et influents
x11()
par(mfrow = c(1,2))

# Points leviers (Hat values) évalués par rapport au seuil 2*p/n
leviers_log = hatvalues(log_reg6)
plot(leviers_log, type="h", main="Points leviers (log_reg6)", xlab="Index", ylab="Hat values")
abline(h = 2 * length(log_reg6$coefficients) / nrow(data_log), col="red")
grid()

# Points influents (Distance de Cook) évalués par rapport au quantile de Fisher
dist_log = cooks.distance(log_reg6)
plot(dist_log, type="h", main="Distance de Cook (log_reg6)", xlab="Index", ylab="Distance")
abline(h = qf(0.1, 3, 17), col="red")
grid()
par(mfrow = c(1,1))

# Les points 4 (2009) et 15 (2020) dépassent les seuils, on les enlève
data_log2 = data_log[-15, ]
data_log2 = data_log2[-4, ]


# Deuxième recherche : modèle sans les années 2009 et 2020
# Élimination pas à pas des variables qui ont leur p-value > 0.05
log_reg21 = lm(TotalPassengers_US~ . - Annee, data_log2)
summary(log_reg21)

# On enlève Investissements (plus grande p-value)
log_reg22 = lm(TotalPassengers_US~ . - Annee - Investissements, data_log2)
summary(log_reg22)

# On enlève CPI
log_reg23 = lm(TotalPassengers_US~ . - Annee - Investissements - CPI, data_log2)
summary(log_reg23)

# On enlève GDP
log_reg24 = lm(TotalPassengers_US~ . - Annee - Investissements - CPI - GDP, data_log2)
summary(log_reg24)

# On enlève Conso_gouv
log_reg25 = lm(TotalPassengers_US~ . - Annee - Investissements - CPI - GDP - Conso_gouv, data_log2)
summary(log_reg25)

# On enlève Exports
log_reg26 = lm(TotalPassengers_US~ . - Annee - Investissements - CPI - GDP - Conso_gouv - Exports, data_log2)
summary(log_reg26)

# On enlève Conso_Privee pour obtenir le modèle réduit
log_reg27 = lm(TotalPassengers_US~ . - Annee - Investissements - CPI - GDP - Conso_gouv - Exports - Conso_Privee, data_log2)
summary(log_reg27)

# Test de Fisher pour valider la réduction du modèle
anova(log_reg21, log_reg27)
# La p-value > 0.05 valide notre modèle réduit

# Analyse des résidus
residu_log2 = log_reg27$residuals
sum(residu_log2)
# La somme des résidus est bien égale à 0

# Test d'autocorrélation de Durbin-Watson
dwtest(log_reg27)
# On a bien p-value > 0.05 donc on ne rejette pas H0

# Affichage des graphiques de diagnostic de base
x11()
par(mfrow = c(2,2))
plot(log_reg27)
par(mfrow = c(1,1))

# Diagnostic des points leviers et influents
x11()
par(mfrow = c(1,2))

# Points leviers (Hat values) évalués par rapport au seuil 2*p/n
leviers_log2 = hatvalues(log_reg27)
plot(leviers_log2, type="h", main="Points leviers (log_reg27)", xlab="Index", ylab="Hat values")
abline(h = 2 * length(log_reg27$coefficients) / nrow(data_log2), col="red")
grid()

# Points influents (Distance de Cook) évalués par rapport au quantile de Fisher
dist_log2 = cooks.distance(log_reg27)
plot(dist_log2, type="h", main="Distance de Cook (log_reg27)", xlab="Index", ylab="Distance")
abline(h = qf(0.1, 2, nrow(data_log2) - 2), col="red")
grid()
par(mfrow = c(1,1))

# L'année 2021 (point 14 après décalage) dépasse le seuil, on l'enlève
data_log3 = data_log2[-14,]


# Troisième recherche : modèle final sans 2009, 2020 et 2021
# Élimination pas à pas des variables qui ont leur p-value > 0.05
log_reg31 = lm(TotalPassengers_US~ . - Annee, data_log3)
summary(log_reg31)

# On enlève Investissements (plus grande p-value)
log_reg32 = lm(TotalPassengers_US~ . - Annee - Investissements, data_log3)
summary(log_reg32)

# On enlève Imports
log_reg33 = lm(TotalPassengers_US~ . - Annee - Investissements - Imports, data_log3)
summary(log_reg33)

# On enlève Conso_Privee pour obtenir le modèle final
log_reg34 = lm(TotalPassengers_US~ . - Annee - Investissements - Imports - Conso_Privee, data_log3)
summary(log_reg34)

# Test de Fisher pour valider la réduction du modèle
anova(log_reg31, log_reg34)
# La p-value > 0.05 valide notre modèle final

# Analyse des résidus
residu_log3 = log_reg34$residuals
sum(residu_log3)
# La somme des résidus est bien égale à 0

# Test d'autocorrélation de Durbin-Watson
dwtest(log_reg34)
# On a bien p-value > 0.05 donc on ne rejette pas H0

# Affichage des graphiques de diagnostic de base
x11()
par(mfrow = c(2,2))
plot(log_reg34)
par(mfrow = c(1,1))

# Validation finale des points leviers et influents
x11()
par(mfrow = c(1,2))

# Points leviers (Hat values) évalués par rapport au seuil 2*p/n
leviers_log3 = hatvalues(log_reg34)
plot(leviers_log3, type="h", main="Points leviers (log_reg34)", xlab="Index", ylab="Hat values")
abline(h = 2 * length(log_reg34$coefficients) / nrow(data_log3), col="red")
grid()

# Points influents (Distance de Cook) évalués par rapport au quantile de Fisher
dist_log3 = cooks.distance(log_reg34)
plot(dist_log3, type="h", main="Distance de Cook (log_reg34)", xlab="Index", ylab="Distance")
abline(h = qf(0.1, 2, nrow(data_log3) - 2), col="red")
grid()
par(mfrow = c(1,1))

# Les valeurs sont proches ou sous les seuils, on les conserve

# Le R2 ajusté (0.9753) confirme la qualité du modèle log_reg34
summary(log_reg34)

# Récupération des coefficients Beta du modèle
betas = log_reg34$coefficients
betas

# Intervalles de confiance des coefficients du modèle final
ic_betas = confint(log_reg34, level = 0.95)
ic_betas

# Tableau final des coefficients betas du modèle avec un IC à 95%
coef_betas_avec_IC = data.frame(
  Borne_Basse = ic_betas[, 1],
  Beta = betas,
  Borne_Haute = ic_betas[, 2]
)

# Affichage du tableau final
coef_betas_avec_IC

# Prévisions pour la période 2026 - 2034
# Chargement et passage en log des données économiques futures
data_futur = data_base[21:29, ]
log_data_futur = data_futur
log_data_futur[, -1] = log(data_futur[, -1])

# Prédiction des passagers en log avec intervalle de confiance à 95%
predictions_log = predict(log_reg34, newdata = log_data_futur, interval = "prediction", level = 0.95)

# Retour aux valeurs réelles avec l'exponentielle
predictions_exp = exp(predictions_log)

# Extraction des données pour le tracé
annees_futures = data_futur$Annee
prevision_centrale = predictions_exp[, "fit"]
borne_basse = predictions_exp[, "lwr"]
borne_haute = predictions_exp[, "upr"]

# Tracé du graphique historique et prévisions
annees_totales = c(data$Annee, annees_futures)
passagers_historiques = data$TotalPassengers_US
passagers_totaux = c(passagers_historiques, prevision_centrale)

plot(annees_totales, passagers_totaux, 
     type = "b",                 
     pch = 19,                  
     col = "blue",              
     lwd = 2,  
     ylim = c(min(passagers_historiques, borne_basse), max(passagers_totaux, borne_haute)),
     xlab = "Années", 
     ylab = "Nombre total de passagers US",
     main = "Trafic aérien : Historique complet et Prévisions (2006 - 2034)")

# Ajout de la ligne de séparation passé/futur
abline(v = 2025.5, col = "darkgrey", lty = 3, lwd = 2)
text(2025.5, max(borne_haute), "Début des prévisions", pos=2, col="darkgrey", cex=0.8)

# Ajout des bornes de l'intervalle de confiance
lines(annees_futures, borne_basse, col = "red", lty = 2, lwd = 1.5)
lines(annees_futures, borne_haute, col = "red", lty = 2, lwd = 1.5)

# Ajout de la légende
legend("topleft", 
       legend = c("Trafic Réel / Prévision", "Intervalle de Confiance (95%)"),
       col = c("blue", "red"), 
       lty = c(1, 2), 
       pch = c(19, NA),
       bty = "n")

