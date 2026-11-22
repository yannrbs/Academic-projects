# Code R1

# Simulation d'un processus de poisson homogène

simul_HPP = function(t, lambda) {
  # Génération du nombre total de sauts sur [0, t], Nt suit une loi de poisson de paramètre lambda * t
  n = rpois(1, lambda * t)
  
  if (n == 0) {
    # Aucun saut
    return(numeric(0)) 
  }
  
  # Génération de n variables uniformes sur [0, t]
  U = runif(n, min = 0, max = t)
  
  # Les temps de sauts sont croissant, on les trient donc du plus petit au plus grand
  jump_times = sort(U)
  
  return(jump_times)
}

set.seed(3)

# Paramètres du modèle
lambda = 2

# Tracé de la trajectoire de Nt 
jump_times_traj = simul_HPP(100, lambda)

# Ajout de l'origine (0,0)
times = c(0, jump_times_traj)
counts = 0:length(jump_times_traj)

par(mfrow = c(1, 1))
plot(times, counts, type = "s", col = "blue", lwd = 2,
     xlab = "Temps (t)", ylab = expression(N[t]),
     main = paste("Trajectoire d'un Processus de Poisson (lambda =", lambda, ")"))
grid()

# Verification que les Sn suivent une loi expo de param lambda
t_valeurs = c(50, 100, 500, 1000)

# Affichage de 4 graphiques 
par(mfrow = c(2, 2))

# Boucle poour tracé les résultats du modèle en fonction de plusieurs t 
for (t in t_valeurs) {
  # Simulation pour le t courant
  jump_times = simul_HPP(t, lambda)
  
  # Calcul des temps de séjour, Sn = Tn - Tn-1
  tps_sejour = diff(c(0, jump_times))
  
  # Récupération du nombre réel de sauts générés
  n_sauts = length(jump_times)
  
  # Tracé de l'histogramme 
  hist(tps_sejour, breaks = 20, prob = TRUE, col = "lightblue", border = "white",
       xlab = "Temps de séjour (Sn)", ylab = "Densité",
       main = paste("t =", t, "|", n_sauts, "sauts"))
  
  # Superposition de la courbe théorique
  curve(dexp(x, rate = lambda), col = "red", lwd = 2, add = TRUE)
  
  # Légende
  legend("topright", legend = paste("Densité théorique : Exp(", lambda, ")", sep=""), 
         col = "red", lwd = 2, bty = "n", cex = 0.9)
  grid()
}

# Plus on augmente le temps t, plus on observe de sauts. 
# Grâce à la Loi Forte des Grands Nombres (LFGN), cet échantillon plus grand 
# permet à l'histogramme de converger vers la courbe théorique d'une loi exp de paramètre lambda
# Cela confirme bien la théorie : dans un processus de Poisson (HPP), les temps de séjour Sn suivent une loi exponentielle de paramètre lambda

# Réinitialisation de la fenêtre graphique
par(mfrow = c(1, 1))


# Code R2

# Simulation d'une trajectoire du système M/M/K sur [0, t]
simuler_mmk = function(x0, K, lambda, mu, t) {
  temps_actuel = 0
  etat_actuel = x0
  
  # Historique des instants de saut et des états occupés
  temps = c(0)
  etats = c(x0)
  
  while (temps_actuel < t) {
    # Nombre de serveurs actifs : au maximum K
    serveurs_actifs = min(etat_actuel, K)
    taux_total = lambda + serveurs_actifs * mu
    
    # Temps de séjour dans l'état actuel (loi exponentielle)
    dt = rexp(1, rate = taux_total)
    
    # Si le saut suivant dépasse l'horizon t, on s'arrête
    if (temps_actuel + dt > t) {
      break
    }
    
    # Mise à jour du temps
    temps_actuel = temps_actuel + dt
    
    # Nature du saut : arrivée (+1) ou fin de service (-1)
    prob_arrivee = lambda / taux_total
    if (runif(1) < prob_arrivee) {
      etat_actuel = etat_actuel + 1
    } else {
      etat_actuel = etat_actuel - 1
    }
    
    temps = c(temps, temps_actuel)
    etats = c(etats, etat_actuel)
  }
  
  # Ajout du point final en t pour fermer proprement la trajectoire sur [0, t]
  temps = c(temps, t)
  etats = c(etats, etat_actuel)
  
  return(data.frame(temps = temps, etats = etats))
}

# Aperçu des 10 premiers sauts d'une trajectoire
traj_exemple = simuler_mmk(x0 = 0, K = 2, lambda = 2, mu = 2, t = 200)
head(traj_exemple, 10)




# Question 2 
# Fonction d'affichage d'une trajectoire
tracer_mmk = function(trajectoire, titre = "Trajectoire du système M/M/K") {
  plot(trajectoire$temps, trajectoire$etats, type = "s",
       col = "darkblue", lwd = 2,
       xlab = "Temps t", ylab = "Nombre de clients dans le système",
       main = titre, ylim = c(0, max(trajectoire$etats) + 1))
  grid()
}

# Illustration des propriétés du modèle selon les paramètres K, lambda, mu et x0
set.seed(20)
t_max = 30

traj_1 = simuler_mmk(x0 = 50,  K = 5, lambda = 4, mu = 2, t = t_max) # Stable (lambda < K*mu)
tracer_mmk(traj_1)

traj_2 = simuler_mmk(x0 = 50,  K = 5, lambda = 5, mu = 1, t = t_max) # critique (lambda = K*mu)
tracer_mmk(traj_2)

traj_3 = simuler_mmk(x0 = 50,  K = 5, lambda = 7, mu = 1, t = t_max) # explosif (lambda > K*mu)
tracer_mmk(traj_3)




#Question 3 


#Estimation de la proba invariante avec la trajectoire
estimer_pi_empirique = function(trajectoire, t) {
  etat_max = max(trajectoire$etats)
  durees = diff(trajectoire$temps)
  etats_visites = trajectoire$etats[-nrow(trajectoire)]
  
  pi_emp = numeric(etat_max + 1)
  for (i in 0:etat_max) {
    pi_emp[i + 1] = sum(durees[etats_visites == i]) / t
  }
  return(pi_emp)
}


#Calcul de la proba invariante théorique
calculer_pi_theorique = function(K, lambda, mu, etat_max) {
  rho = lambda / (K * mu)
  if (rho >= 1) stop("Condition de stabilité non respectée : lambda doit être < K * mu")
  
  a = lambda / mu
  k = 0:(K - 1)
  pi_0 = 1 / (sum((a^k) / factorial(k)) + ((a^K) / factorial(K)) * (1 / (1 - rho)))
  
  pi_theo = numeric(etat_max + 1)
  for (i in 0:etat_max) {
    if (i < K) {
      pi_theo[i + 1] = pi_0 * (a^i) / factorial(i)
    } else {
      pi_theo[i + 1] = pi_0 * (a^K) * (rho^(i - K)) / factorial(K)
    }
  }
  return(pi_theo)
}



# Comparaison de la proba inv théorique et empirique

comparer_pi_mmk = function(x0, K, lambda, mu, t) {
  # Simulation de la trajectoire avec la fonction simuler_mmk
  traj = simuler_mmk(x0, K, lambda, mu, t)
  etat_max = max(traj$etats)
  
  #Calcul des distributions empirique et théorique
  pi_emp  = estimer_pi_empirique(traj, t)
  pi_theo = calculer_pi_theorique(K, lambda, mu, etat_max)
  
  # Matrice pour le barplot (affichage limité aux états <= 20 pour la lisibilité)
  limite_affichage = min(etat_max, 20)
  mat_comparaison = rbind(
    Empirique = pi_emp[1:(limite_affichage + 1)],
    Theorique = pi_theo[1:(limite_affichage + 1)]
  )
  colnames(mat_comparaison) = 0:limite_affichage
  
  # Tracé du diagramme en barres juxtaposées
  barplot(mat_comparaison, beside = TRUE,
          col = c("blue", "red"),
          xlab = "Nombre de clients dans le système (état i)",
          ylab = "Probabilité invariante",
          main = paste("Distribution empirique vs Probabilité invariante (t =", t, ")"),
          legend.text = c("Estimation empirique", "Probabilité invariante théorique (pi)"),
          args.legend = list(x = "topright", bg = "white"))
  grid(nx = NA, ny = NULL) # Grille horizontale uniquement
  
  # Évaluation de la convergence (Erreur en norme L1)
  erreur_L1 = sum(abs(pi_emp - pi_theo))
  cat("Erreur en norme L1 ||pi_emp - pi_theo||_1 :", round(erreur_L1, 5), "\n")
  
  return(invisible(erreur_L1))
}

# Executions
set.seed(20)

#Pour t_max = 5
comparer_pi_mmk(x0 = 20, K = 5, lambda = 4, mu = 2, t = 5)
#Pour t_max = 5000
comparer_pi_mmk(x0 = 20, K = 5, lambda = 4, mu = 2, t = 5000)





# Validation rigoureuse par test du Khi-deux 

set.seed(2)

# Choix aléatoire des paramètres en garantissant la stabilité du système (lambda < K * mu)
K = sample(1:5, 1)
mu = runif(1, min = 1, max = 4)
lambda = runif(1, min = 1, max = K * mu - 0.1) 
t_long = 5000 # Horizon temporel suffisamment grand

# Simulation et récupération des distributions
traj_alea = simuler_mmk(x0 = 0, K = K, lambda = lambda, mu = mu, t = t_long)
etat_max = max(traj_alea$etats)

pi_emp = estimer_pi_empirique(traj_alea, t_long)
pi_theo = calculer_pi_theorique(K, lambda, mu, etat_max)

# Test d'adéquation du Khi-deux
# On convertit les temps de séjour en pseudo-comptages via floor(pi_emp * t_long)
comptages_observes = floor(pi_emp * t_long)

# Normalisation de la loi theorique tronque
pi_tronque = pi_theo / sum(pi_theo)

# Utilisation de la fonction chi-deux avec p-valeur simulée sur 5 000 itérations
res_test = chisq.test(x = comptages_observes, p = pi_tronque, 
                      simulate.p.value = TRUE, B = 5000)
print(res_test)



# Code R3

# Question 1

# Fonction pour calculer pi_t et retourner aussi le nombre d'itérations k
state_vector_pi_t = function(t, lambda_hat, Q, alpha, error){
  # Définition de P_hat
  P_hat = diag(nrow(Q)) + (1/lambda_hat) * Q
  
  # Initialisation de pi_t
  pi_t = rep(0, length(alpha))
  
  # Initialisation de alpha_k
  alpha_k = alpha
  
  # Initialisation de k 
  k = 0
  
  # Paramètre de la loi de Poisson
  mu = lambda_hat * t
  
  repeat {
    # Calcul de pi_t
    pi_t = pi_t + dpois(k, mu) * alpha_k
    
    # Actualisation de alpha * P_hat^k
    alpha_k = alpha_k %*% P_hat
    
    # Calcul de la borne d'erreur Rn théorique démontrée à la Q3 sur l'uniformisation des JMP finis
    # Y suit une loi de Poisson, ppois(k, mu) = P(Y <= k), donc 1 - ppois = P(Y > k) 
    Rn = 1 - ppois(k, mu)
    
    if (Rn < error){
      break
    }
    
    # Actualisation de k 
    k = k + 1
  }
  
  # On retourne une liste contenant le vecteur pi_t et le nombre d'itérations k
  return(list(pi_t = pi_t, iterations = k))
}

set.seed(3)


# Question 2

# Test 1 sur un générateur simple

# Définition d'un générateur Q simple (2 états)
Q1 = matrix(c(-2,  2,  
              3, -3), nrow = 2, byrow = TRUE)

# Paramètres du test 1
# lambda_hat1 doit être strictement supérieur à max(|Q1(i,i)|), ici 3
lambda_hat1 = 4  

# On démarre avec certitude dans l'état 1
alpha1 = c(1, 0) 

# Un temps assez grand pour converger
t1 = 10      

# Précision de l'erreur
error1 = 1e-6   

# Exécution de la fonction pour le test 1
res_test1 = state_vector_pi_t(t1, lambda_hat1, Q1, alpha1, error1)
pi_final1 = res_test1$pi_t

# Affichage test 1
cat("Vecteur d'état estimé pi_t (Test 1) :", pi_final1, "\n")
cat("Probabilité invariante théorique : 0.6 0.4 \n\n")

# Détermination du temps t* (fin du régime transitoire)
# Séquence de temps (de 0 à 10, par pas de 0.01 pour une bonne précision sur t*)
t_seq_test1 = seq(0, 10, by = 0.01)

# Matrice pour stocker les résultats (chaque ligne = un instant t, chaque colonne = un état)
pi_matrix_test1 = matrix(0, nrow = length(t_seq_test1), ncol = 2)

# Calcul de pi_t pour chaque instant t de la séquence
for (i in 1:length(t_seq_test1)) {
  pi_matrix_test1[i, ] = state_vector_pi_t(t_seq_test1[i], lambda_hat1, Q1, alpha1, error1)$pi_t
}

# pi théorique pour Q1
pi_theo1 = c(0.6, 0.4)

# Choix d'un seuil global eta basé sur une erreur relative de 0.5% sur l'état le plus petit
eta1 = 0.005 * min(pi_theo1) 

# Calcul de la distance infinie max_i |pi_t(i) - pi_theo1(i)| à chaque instant t
ecart_max1 = apply(abs(pi_matrix_test1 - matrix(pi_theo1, nrow = length(t_seq_test1), ncol = 2, byrow = TRUE)), 1, max)

# Détection de l'instant t* où la distance infinie passe sous eta1
indice_conv1 = which(ecart_max1 <= eta1)[1]
t_star1 = t_seq_test1[indice_conv1]

cat("Seuil de tolérance stationnaire (eta1) :", eta1, "\n")
cat("Temps d'entrée en régime stationnaire (t*) : t* =", t_star1, "\n")

# Comparaison de l'efficacité de différents lambda_hat
# On définit un vecteur de lambda allant de 3.1 à 10 par pas de 0.1
lambda_values = seq(3.1, 10, 0.1)

# On définit un dataframe pour stocker les résultats
resultats_comparaison = data.frame(lambda_hat = numeric(), iterations = numeric(), pi_t_1 = numeric())

# On teste chaque valeur de lambda_hat dans notre fonction state_vector_pi_t
for (lam in lambda_values) {
  res = state_vector_pi_t(t1, lam, Q1, alpha1, error1)
  resultats_comparaison = rbind(resultats_comparaison, data.frame(lambda_hat = lam, iterations = res$iterations, pi_t_1 = res$pi_t[1]))
}

# Tri des résultats par nombre d'itérations k croissant
resultats_tries = resultats_comparaison[order(resultats_comparaison$iterations), ]

# Graphique du nombre d'itérations k en fonction des différentes valeurs de lambda_hat
plot(resultats_tries$lambda_hat, resultats_tries$iterations, 
     type = "l", lwd = 2, col = "blue",
     xlab = "Valeurs de lambda_hat", 
     ylab = "Nombre d'itérations (k)",
     main = "Impact de lambda_hat sur le nombre d'itérations")
grid()

# On voit qu'on doit choisir un lambda_hat le plus proche de max(|Q1(i,i)|)

# Test 2 sur un processus de Naissance et de Mort à 10 états

# Construction du générateur Q2 (Processus de Naissance et de Mort à N2 états)
N2 = 10
Q2 = matrix(0, nrow = N2, ncol = N2)

# On choisit les taux de naissance et de mort
taux_naissance2 = 2
taux_mort2 = 2

for (i in 1:N2) {
  if (i > 1) Q2[i, i-1] = taux_mort2        # Diagonale inférieure (Mort)
  if (i < N2) Q2[i, i+1] = taux_naissance2 # Diagonale supérieure (Naissance)
  Q2[i, i] = -sum(Q2[i, ])                 # La somme de la ligne doit faire 0
}

# Paramètres du test 2
# Le terme diagonal maximum en valeur absolue est |-4|
# On choisit donc un lambda_hat2 strictement supérieur et proche de 4, ici 4.1
lambda_hat2 = 4.1  

# On démarre avec proba 1 dans le premier état
alpha2 = c(1, rep(0, N2 - 1)) 

# Temps d'observation assez long pour converger vers la loi invariante
t2 = 60         
error2 = 1e-6    

# Exécution de la fonction pour le test 2 (extraction de $pi_t)
res_test2 = state_vector_pi_t(t2, lambda_hat2, Q2, alpha2, error2)
pi_final2 = res_test2$pi_t

# Affichage test 2
cat("Vecteur d'état estimé pi_t (Test 2) :\n")
print(pi_final2)

cat("\nProbabilité invariante théorique (Loi uniforme) :\n")
print(rep(0.1, N2))


# Illustration de la convergence de pi_t vers la loi invariante et le changement de régime du processus
# (Basée sur les paramètres du Test 2)

# Séquence de temps (de 0 à 60, par pas de 0.01 pour une bonne précision sur t*)
t_seq2 = seq(0, 60, by = 0.01)

# Matrice pour stocker les résultats (chaque ligne = un instant t, chaque colonne = un état)
pi_t_matrix2 = matrix(0, nrow = length(t_seq2), ncol = N2)

# Calcul de pi_t pour chaque instant t de la séquence (extraction de $pi_t)
for (i in 1:length(t_seq2)) {
  pi_t_matrix2[i, ] = state_vector_pi_t(t_seq2[i], lambda_hat2, Q2, alpha2, error2)$pi_t
}

# Détermination du seuil eta et du temps de convergence t*
pi_theo2 = 0.1

# Choix de eta : erreur relative de 0.5% autour de 0.1 (et bien supérieur à error2 = 1e-6)
eta = pi_theo2 * 0.005

# Calcul de la distance infinie max_i |pi_t(i) - 0.1| à chaque instant t
ecart_max2 = apply(abs(pi_t_matrix2 - pi_theo2), 1, max)

# Détection de l'instant t* où toutes les courbes entrent dans le couloir [0.1 - eta, 0.1 + eta]
indice_conv = which(ecart_max2 <= eta)[1]
t_star = t_seq2[indice_conv]

cat("Seuil de tolérance stationnaire (eta) :", eta, "\n")
cat("Fin du régime transitoire (entrée dans le couloir +/- eta) à t* =", t_star, "\n")


# Graphique 1 : Vue globale (Régime transitoire puis stationnaire) 
matplot(t_seq2, pi_t_matrix2, type = "l", lty = 1, lwd = 2, col = 1:N2,
        xlab = "Temps t", ylab = "Probabilité P(X_t = i)",
        main = "Convergence du vecteur d'état vers la probabilité invariante",
        ylim = c(0, 1))

# Ligne horizontale de la limite théorique (0.1)
abline(h = pi_theo2, col = "black", lty = 2, lwd = 1.5)

# Ligne verticale marquant la frontière entre régime transitoire (t < t*) et stationnaire (t >= t*)
abline(v = t_star, col = "red", lty = 3, lwd = 2)
text(t_star + 1, 0.8, paste("Fin du régime transitoire\nt* =", t_star), col = "red", adj = 0, cex = 0.85)

legend("topright", legend = paste("État", 1:N2), col = 1:N2, lty = 1, cex = 0.7, ncol = 2)
grid()

# Graphique 2 : Zoom sur l'entrée dans le couloir +/- eta (t de 25 à 60) 
indices_zoom = which(t_seq2 >= 25)

matplot(t_seq2[indices_zoom], pi_t_matrix2[indices_zoom, ], type = "l", lty = 1, lwd = 2, col = 1:N2,
        xlab = "Temps t", ylab = "Probabilité P(X_t = i)",
        main = paste("Zoom sur l'entrée en régime stationnaire (couloir +/-", eta, ")"),
        ylim = c(0.0985, 0.1015)) # Fenêtre centrée autour du couloir [0.0995, 0.1005]

# Ligne théorique centrale (0.1)
abline(h = pi_theo2, col = "black", lty = 1, lwd = 1.5)

# Tracé du couloir de convergence [0.1 - eta, 0.1 + eta] en pointillés rouges
abline(h = pi_theo2 + eta, col = "red", lty = 2, lwd = 2)
abline(h = pi_theo2 - eta, col = "red", lty = 2, lwd = 2)

# Ligne verticale bleue indiquant l'instant exact t* où la dernière courbe entre dans le couloir
abline(v = t_star, col = "blue", lty = 3, lwd = 2)
text(t_star + 0.5, 0.1012, paste("t* =", t_star), col = "blue", adj = 0, font = 2)

legend("topright", 
       legend = c(paste("État", 1:N2), "Couloir +/- eta", "Temps conv. t*"),
       col = c(1:N2, "red", "blue"), lty = c(rep(1, N2), 2, 3), lwd = 1.5, cex = 0.65, ncol = 2)
grid()



# Code R4

# Question 1

# Partie S1, loi du nombre d'évenements
# On fait une fonction pour simuler n evenements

simul_n_event = function(a, b, t, lambda_hat, Q, error){
  # Définition de P_hat
  P_hat = diag(nrow(Q)) + (1/lambda_hat) * Q
  
  # Nombre d'états
  N_states = nrow(Q)
  
  # Definition de alpha_a 
  alpha_a = rep(0, N_states)
  
  # On part de l'état a avec proba 1
  alpha_a[a] = 1
  
  # Calcul du dénominateur exp(Qt)(a, b) 
  # Revient à calculer P(Xt = b | X0 = a) = pi_t avec X0 ~ a
  # On prend donc l'element b de notre fonction state_vector_pi_t
  denom = state_vector_pi_t(t, lambda_hat, Q, alpha_a, error)$pi_t[b]
  
  # Cas ou on ne peut pas atteindre l'état b depuis l'état a
  if (denom <= 0){
    stop("La probabilité d'atteindre l'état b depuis l'état a est nulle.")
  }
  
  # Cas ou c'est possible
  # Calcul du numérateur avec la méthode d'inversion 
  
  # On tire U sur [0, 1]
  u = runif(1, 0, 1)
  
  # On initialise la probabilité cumulée F(n)
  cum_prob = 0 
  
  # On commence à n = 0
  n = 0  
  
  # On initialise P_hat^n, qui vaut l'identité car P_hat^0 = Id
  P_hat_n = diag(N_states)
  
  # On crée une liste pour stocker les P_hat^n
  P_powers = list()
  
  # Calcul du numérateur
  repeat{
    # Stockage de la matrice P_hat^n actuel
    P_powers[[n+1]] = P_hat_n
    
    # Calcul du numerateur 
    num = P_hat_n[a, b] * dpois(n, lambda_hat * t) / denom
    
    cum_prob = cum_prob + num
    
    # Si u <= cum_prob alors on sort de la boucl et on renvoie n le nombre d'évenements 
    if (u <= cum_prob){
      break
    }
    
    # On met à jour P_hat^n
    P_hat_n = P_hat_n %*% P_hat
    
    # Mise à jour de n 
    n = n + 1
  }
  
  return(list(n, P_powers))
}


# Fonction pour simuler une JMP qui commence en a et fini en b 

simul_JMP_ab = function(a, b, t, lambda_hat, Q, error){
  N_states = nrow(Q)
  
  # Simulation de n et récupération des puissances de P_hat (S1)
  res_s1 = simul_n_event(a, b, t, lambda_hat, Q, error)
  n = res_s1[[1]]
  P_powers = res_s1[[2]]
  
  # Cas ou n = 0 (S2)
  if (n == 0) {
    # Le processus reste dans l'état a de 0 jusqu'à t
    return(list(times = c(0, t), states = c(a, a)))
  }
  
  # Cas ou n = 1 (S3)
  if (n == 1) {
    if (a == b) {
      # Cas a = b : pas de changement d'état effectif
      return(list(times = c(0, t), states = c(a, a)))
    } else {
      # Cas a != b : un unique saut à un instant t1 tiré uniformément sur [0, t]
      t1 = runif(1, 0, t)
      return(list(times = c(0, t1, t), states = c(a, b, b)))
    }
  }
  
  # Cas n >= 0, génération des n instants d'événements triés (S4)
  # On fait n tirage uniforme entre 0 et t avec n qui provient de simul_n_event (S4.1)
  t_events = sort(runif(n, 0, t))
  
  # Simulation des états intermédiaires z_1, ..., z_{n-1} via la formule (4) (S4.2 et S4.3)
  z = rep(0, n + 1)
  
  # Z_0 = a, état initial fixé 
  z[0 + 1] = a 
  
  # Z_n = b, état final fixé
  z[n + 1] = b
  
  # On calcul P_hat
  P_hat = diag(N_states) + (1 / lambda_hat) * Q
  
  # Calcul de la formule (4)
  if (n >= 2) {
    for (i in 1:(n - 1)) {
      # On récupère l'état i - 1, z_{i-1}
      z_prev = z[i - 1 + 1]
      
      # Calcul de la loi conditionnelle de Z_i selon la formule (4) du sujet
      # P(Z_i = z | Z_{i-1} = z_prev, Z_n = b)
      
      # P_hat(z_prev, z) * P_hat^{n-i}(z, b)
      num_vec = P_hat[z_prev, ] * P_powers[[n - i + 1]][, b] 
      
      # P_hat^{n-i+1}(z_prev, b)
      den_val = P_powers[[n - i + 1 + 1]][z_prev, b]
      
      # Calcul des probabilités conditionnelles pour chaque état possible
      proba_cond = num_vec / den_val
      
      # Sécurité pour ne pas diviser par 0
      proba_cond[is.na(proba_cond)] = 0
      
      # Normalisation
      proba_cond = proba_cond / sum(proba_cond) 
      
      # Tirage de l'état z_i selon la loi que l'on vient de calculer : proba_cond
      z[i + 1] = sample(1:N_states, size = 1, prob = proba_cond)
    }
  }
  
  # Élimination des transitions virtuelles (où z_i = z_{i-1}) (S4.4)
  times_init = c(0, t_events)
  states_init = z 
  
  # On garde l'instant 0 et uniquement les instants où l'état change réellement
  real_jumps = c(TRUE, states_init[2:(n + 1)] != states_init[1:n])
  
  # On rajoute l'instant final t pour fermer le dernier palier sur le graphique
  final_times = c(times_init[real_jumps], t)
  final_states = c(states_init[real_jumps], b)
  
  return(list(times = final_times, states = final_states))
}

# Fonction pour tracer la trajectoire conditionnée 

plot_JMP_bridge = function(res_simul, a, b, t) {
  # Utilisation de type = "s" pour tracer en escalier
  plot(res_simul$times, res_simul$states, type = "s", lwd = 2, col = "blue",
       xlab = "Temps", ylab = "État occupé",
       main = paste("Pont Markovien : X_0 =", a, " et X_t =", b),
       xlim = c(0, t), 
       ylim = c(min(res_simul$states) - 0.5, max(res_simul$states) + 0.5),
       yaxt = "n") # On désactive l'axe y par défaut pour le refaire proprement
  
  # Axe Y avec uniquement des entiers (les états)
  axis(2, at = min(res_simul$states):max(res_simul$states))
  
  # Ajout de gros points rouges pour bien prouver que les conditions aux bornes sont respectées
  points(0, a, col = "red", pch = 16, cex = 1.5)
  points(t, b, col = "red", pch = 16, cex = 1.5)
  grid()
}

# Test 1 : Sur le générateur Q1 (2 états) 
# On veut simuler un chemin qui part de 1 et finit en 2 au temps t = 5
set.seed(42)

# Paramètres 
a1 = 1
b1 = 2
t1_bridge = 5
lambda_hat1 = 3.1
error = 1e-6

# Simulation et affichage
res_bridge1 = simul_JMP_ab(a = a1, b = b1, t = t1_bridge, lambda_hat = lambda_hat1, Q = Q1, error = error)
plot_JMP_bridge(res_bridge1, a = a1, b = b1, t = t1_bridge)


# Test 2 : Sur le générateur Q2 (10 états, Naissance et Mort)
# On part de l'état 1 et on force le processus à être dans l'état 8 au temps t = 10
set.seed(90) 

# Paramètres
a2 = 1
b2 = 8
t2_bridge = 10
lambda_hat2 = 4.1

# Simulation et affichage
res_bridge2 = simul_JMP_ab(a = a2, b = b2, t = t2_bridge, lambda_hat = lambda_hat2, Q = Q2, error = error)
plot_JMP_bridge(res_bridge2, a = a2, b = b2, t = t2_bridge)


# Test 3 : Cas où a = b sur Q2 
# Le processus part de 3 et finit en 3 (il doit faire des allers-retours ou rester sur place)
set.seed(9)

# Paramètres
a3 = 3
b3 = 3
t3_bridge = 10

# Simulation et affichage
res_bridge3 = simul_JMP_ab(a = a3, b = b3, t = t3_bridge, lambda_hat = lambda_hat2, Q = Q2, error = error)
plot_JMP_bridge(res_bridge3, a = a3, b = b3, t = t3_bridge)



# Question 2 

# Fonction pour simuler une JMP avec des points intermediaires fixes
simul_JMP_points_fixes = function(a, b, t, t_points, x_points, lambda_hat, Q, error) {
  # Création des vecteurs complets incluant le départ et l'arrivée
  T_full = c(0, t_points, t)
  X_full = c(a, x_points, b)
  
  # Sécurité, on vérifie que les tailles correspondent
  if (length(t_points) != length(x_points)) {
    stop("Le nombre de temps intermédiaires doit être égal au nombre d'états intermédiaires.")
  }
  
  # Initialisation des vecteurs finaux
  times_final = c()
  states_final = c()
  
  # Boucle sur chaque sous-intervalle
  for (i in 1:(length(T_full) - 1)) {
    # Extraction des paramètres du sous-intervalle
    t_start = T_full[i]
    t_end = T_full[i + 1]
    
    state_start = X_full[i]
    state_end = X_full[i + 1]
    
    # Durée du sous-intervalle
    delta_t = t_end - t_start
    
    # Simulation du pont markovien sur cette durée (Question 1)
    res_sub = simul_JMP_ab(state_start, state_end, delta_t, lambda_hat, Q, error)
    
    # Décalage des temps simulés pour les replacer dans le bon chronométrage global
    sub_times = res_sub$times + t_start
    sub_states = res_sub$states
    
    # Suppression du premier point pour éviter les doublons aux raccords
    # (sauf pour le tout premier intervalle qui doit inclure le temps 0)
    if (i > 1) {
      sub_times = sub_times[-1]
      sub_states = sub_states[-1]
    }
    
    # Ajout du tronçon à la trajectoire globale
    times_final = c(times_final, sub_times)
    states_final = c(states_final, sub_states)
  }
  
  return(list(times = times_final, states = states_final))
}


# Test 1 : Quelques points de passage
# On utilise le générateur Q2 (10 états) et nos paramètres précédents
set.seed(123)
error = 1e-6
lambda_hat2 = 4.1

# Paramètres du chemin global
a_multi = 1           # Départ en 1
b_multi = 10          # Arrivée en 10
t_fin = 15            # Temps total

# Points de passage obligatoires
# À t=5 on veut être en état 4, à t=10 on veut être en état 8 et à t=12 on veut être en 3
t_inter = c(5, 10, 12)
x_inter = c(4, 8, 3)

# Lancement de la simulation
res_multi = simul_JMP_points_fixes(a_multi, b_multi, t_fin, t_inter, x_inter, lambda_hat2, Q2, error)

# Tracé du graphique complet
plot(res_multi$times, res_multi$states, type = "s", lwd = 2, col = "blue",
     xlab = "Temps t", ylab = "État occupé",
     main = "Trajectoire avec points intermédiaires imposés",
     xlim = c(0, t_fin), 
     ylim = c(1, 10),
     yaxt = "n")

# Refonte de l'axe y avec les états
axis(2, at = 1:10)
grid()

# Ajout des gros points rouges pour prouver que la chaîne passe bien par toutes les contraintes
points(0, a_multi, col = "red", pch = 16, cex = 1.5)               # Départ
points(t_inter, x_inter, col = "red", pch = 16, cex = 1.5)         # Intermédiaires
points(t_fin, b_multi, col = "red", pch = 16, cex = 1.5)           # Arrivée


# Test 2 : Trajectoire contrainte à chaque instant entier
set.seed(404)

# Paramètres du chemin global
a_dense = 1            # Départ en 1 (à t=0)
b_dense = 5            # Arrivée en 5 (à t=20)
t_fin_dense = 20       # Temps total

# On fixe un point de passage pour tous les instants entiers de 1 à 19
t_inter_dense = 1:19

# On choisit une séquence d'états par lesquels le processus doit passer
x_inter_dense = c(2, 3, 5, 7, 9, 10, 8, 6, 4, 2, 1, 3, 4, 6, 8, 9, 7, 6, 4)

# Lancement de la simulation
res_dense = simul_JMP_points_fixes(a_dense, b_dense, t_fin_dense, t_inter_dense, x_inter_dense, lambda_hat2, Q2, error)

# Tracé du graphique complet
plot(res_dense$times, res_dense$states, type = "s", lwd = 2, col = "blue",
     xlab = "Temps t", ylab = "État occupé",
     main = "Trajectoire contrainte à chaque pas de temps entier",
     xlim = c(0, t_fin_dense), 
     ylim = c(1, 10),
     yaxt = "n")

# Refonte de l'axe y avec les états
axis(2, at = 1:10)
grid()

# Ajout des gros points rouges pour prouver que la chaîne passe bien par toutes les contraintes
points(0, a_dense, col = "red", pch = 16, cex = 1.2)
points(t_inter_dense, x_inter_dense, col = "red", pch = 16, cex = 1.2)
points(t_fin_dense, b_dense, col = "red", pch = 16, cex = 1.2)

