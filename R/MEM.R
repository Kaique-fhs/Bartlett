#' @title Bartlett correction in models with errors in variables
#' @name MEM
#'
#' @description Bartlett correction in models with errors in variables
#'
#' @param formula Declare x and y as a vector of numbers of the same size.
#' @param lambda_x Declare as a single number greater than zero.
#' @param lambda_e Declare as a single number greater than zero.
#' @param beta_til Declare as a single number.
#' @param Correction Declare to be TRUE or FALSE. If TRUE, the output is Bartlett corrected.
#'
#'
#' @details The following conditions must be respected:
#'     1 - Declare x and y as a vector of numbers of the same size.
#'     2 - Declare beta_til as a number, remembering that H_0: beta = beta_til.
#'     3 - Declare beta_til as a number, remembering that H_0: beta = beta_til.
#'     4 - Declare lambda_x or lambda_e ratios as a number greater than zero.
#'     5 - Remember to declare only one or the other. Use the one where you will admit to being known.
#'     5 - Default: lambda_x = FALSE, lambda_e = FALSE, beta_til = 0, Correction = TRUE.
#'
#' @return Application of Bartlett's correction in the model with errors in the variables.
#'
#' @author Kaique Souza e Tatiane Silva
#'
#' @seealso \code{\link[stats]{lm}}, \code{\link[base]{+}}
#'
#' @examples
#' x = c(9, 6.6, 12.3, 11.9, 11.9, 12, 9.6, 7.5, 10.9, 10.4, 10.2, 7.4, 11 ,11.8, 8.2)
#' y = c(8, 6, 9.8, 10.8, 9.7, 9.3, 9.2, 6.9, 8.1, 8.7, 8.7, 7.4, 10.1, 10, 7.3)
#' MEM( y~x, lambda_e = 6^(-1), Correction = TRUE)
#'
#' @export
MEM = function(formula = y~x, lambda_x = F, lambda_e = F, beta_til = 1, Correction = T){

  # Condicoes iniciais
  if(is.vector(x) & is.vector(y) & is.numeric(x) & is.numeric(y) & (length(x) == length(y)) &
     (length(lambda_x)==1) & is.logical(Correction) & (length(Correction)==1) &
     (length(lambda_e)==1) & is.numeric(beta_til) & (length(beta_til)==1)){

    ###############################################################
    # Declarando variaveis iniciais
    ### Tamanho da amostra
    n = length(x)

    ##### Inicializacao de vetores, matrizes, arrays e listas
    mu_hat = array(NA,c(2,1))
    mu_til = array(NA,c(2,1))

    Sigma_hat     = array(NA,c(2,2))
    Inv_Sigma_hat = array(NA,c(2,2))

    Sigma_til     = array(NA,c(2,2))
    Inv_Sigma_til = array(NA,c(2,2))

    Z             = array(NA,c(2,1,n))
    diferenca_hat = array(NA,c(2,1,n))
    diferenca_til = array(NA,c(2,1,n))

    Sigma      = list()
    Sigma[[1]] = array(0,c(2,2,5))
    Sigma[[2]] = array(0,c(2,2,5,5))
    Sigma[[3]] = array(0,c(2,2,5,5,5))
    Sigma[[4]] = array(0,c(2,2,5,5,5,5))

    Mu      = list()
    Mu[[1]] = array(0,c(2,1,5))
    Mu[[2]] = array(0,c(2,1,5,5))
    Mu[[3]] = array(0,c(2,1,5,5,5))
    Mu[[4]] = array(0,c(2,1,5,5,5,5))

    D_Inv_Sigma  = list()
    D_Inv_Sigma[[1]] = array(0,c(2,2,5))
    D_Inv_Sigma[[2]] = array(0,c(2,2,5,5))
    D_Inv_Sigma[[3]] = array(0,c(2,2,5,5,5))
    D_Inv_Sigma[[4]] = array(0,c(2,2,5,5,5,5))

    A = array(0,c(5,5,5,5))
    P = array(0,c(5,5,5))
    Q = array(0,c(5,5,5))

    K     = array(NA,c(5,5))
    L     = array(0,c(5,5))
    M     = array(0,c(5,5))
    N     = array(0,c(5,5))
    Inv_K = array(NA,c(5,5))

    K_p_1     = array(NA,c(4,4))
    L_p_1     = array(NA,c(4,4))
    M_p_1     = array(NA,c(4,4))
    N_p_1     = array(NA,c(4,4))
    Inv_K_p_1 = array(NA,c(4,4))

    quantil_90 = stats::qchisq(0.90,1)
    quantil_95 = stats::qchisq(0.95,1)
    quantil_99 = stats::qchisq(0.99,1)

    Produto_hat = matrix(0,nrow=n,ncol=1)
    Produto_til = matrix(0,nrow=n,ncol=1)

    ##### Calculando as somas de quadrados
    x_barra = mean(x)
    y_barra = mean(y)
    S_xx    = mean((x - x_barra)^2)
    S_yy    = mean((y - y_barra)^2)
    S_xy    = mean((x - x_barra)*(y - y_barra))

    ##########################################################################
    ##########################################################################
    ### Fim da declaracao de variaveis

    ###################################################################
    # Comeco do caso lambda_e
    if(lambda_x == F & lambda_e > 0){
      # cat("Lambda_e")

      ##Estimativas de maxima verossimilhanca (e.m.v.)
      beta_hat     = (S_yy - (S_xx*lambda_e) + ( (S_yy - (S_xx*lambda_e))^(2) + 4*lambda_e*(S_xy^(2)) )^(1/2) )/(2*S_xy)
      mu_x_hat     = x_barra
      alfa_hat     = y_barra - (beta_hat*x_barra)
      sigma2_u_hat = ( S_yy - 2*beta_hat*S_xy + beta_hat*beta_hat*S_xx )/( beta_hat*beta_hat + lambda_e )
      sigma2_x_hat = S_xy/beta_hat

      ##Estimativas de maxima verossimilhanca (e.m.v.) sob H0
      alfa_til     = y_barra - (beta_til *x_barra)
      mu_x_til     = ( y_barra*beta_til + lambda_e*x_barra - alfa_til*beta_til )/( lambda_e + beta_til*beta_til)
      sigma2_u_til = ( S_yy - 2*beta_til*S_xy + beta_til*beta_til*S_xx + (x_barra*beta_til - y_barra + alfa_til )^2)/( beta_til*beta_til + lambda_e)
      sigma2_x_til = ( S_yy*(beta_til*beta_til - lambda_e) + 4*beta_til*lambda_e*S_xy - (S_xx*lambda_e*(beta_til*beta_til - lambda_e)) - lambda_e*((x_barra*beta_til + alfa_til)^2) + y_barra*lambda_e*(2*alfa_til -y_barra+(2*x_barra*beta_til)))/( (beta_til*beta_til+lambda_e)^2)

      ##### Matrix de Mdeias
      mu_hat = matrix(c(alfa_hat + (beta_hat*mu_x_hat), mu_x_hat), ncol=1, nrow=2)
      mu_til = matrix(c(alfa_til + (beta_til    *mu_x_til), mu_x_til), ncol=1, nrow=2)

      ##### Matrix de Covariancias
      Sigma_hat = matrix(c((beta_hat*beta_hat*sigma2_x_hat) + (lambda_e*sigma2_u_hat),
                           beta_hat*(sigma2_x_hat),
                           beta_hat*(sigma2_x_hat),
                           sigma2_x_hat + sigma2_u_hat), ncol=2,nrow=2)


      Inv_Sigma_hat = solve(Sigma_hat)
      Det_Sigma_hat = det(Sigma_hat)

      Sigma_til = matrix(c((beta_til*beta_til*sigma2_x_til) + (lambda_e*sigma2_u_til), beta_til*(sigma2_x_til), beta_til*(sigma2_x_til), sigma2_x_til + sigma2_u_til), ncol=2,nrow=2)

      Inv_Sigma_til = solve(Sigma_til)
      Det_Sigma_til = det(Sigma_til)

      ##### Preenchendo os vetores Z, diferenca e produto
      for (i in 1:n)
      {
        Z[,,i]             = matrix(c(y[i],x[i]),ncol = 1,nrow = 2)

        diferenca_hat[,,i] = matrix(Z[,,i] - mu_hat,ncol = 1)
        diferenca_til[,,i] = matrix(Z[,,i] - mu_til,ncol = 1)

        Produto_hat[i]     = (t(matrix(diferenca_hat[,,i],ncol = 1)) %*% Inv_Sigma_hat[,] %*% matrix(diferenca_hat[,,i],ncol = 1))
        Produto_til[i]     = (t(matrix(diferenca_til[,,i],ncol = 1)) %*% Inv_Sigma_til[,] %*% matrix(diferenca_til[,,i],ncol = 1))
      }



      ###################################################################
      # Correction == TRUE - lambda_e
      if(Correction == TRUE){
        # cat("Lambda_e - TRUE")
        ###################################################################
        # Correction == TRUE - lambda_e

        ##### Primeira derivada do vetor mu
        Mu[[1]][1,1,1] = 1
        Mu[[1]][1,1,2] = mu_x_til
        Mu[[1]][1,1,3] = beta_til
        Mu[[1]][2,1,3] = 1


        ##### Segunda derivada do vetor mu
        Mu[[2]][1,1,2,3] = Mu[[2]][1,1,3,2] = 1

        ##### Primeira derivada da matriz Sigma
        Sigma[[1]][1,1,2] = 2*beta_til*sigma2_x_til
        Sigma[[1]][1,2,2] = Sigma[[1]][2,1,2] = sigma2_x_til
        Sigma[[1]][1,1,4] = beta_til^2
        Sigma[[1]][1,2,4] = Sigma[[1]][2,1,4] = beta_til
        Sigma[[1]][2,2,4] = 1
        Sigma[[1]][1,1,5] = lambda_e
        Sigma[[1]][2,2,5] = 1

        ##### Segunda derivada da matriz Sigma
        Sigma[[2]][1,1,2,2]   = 2*sigma2_x_til
        Sigma[[2]][1,1,4,2]   = Sigma[[2]][1,1,2,4] = 2*beta_til
        Sigma[[2]][1,2,4,2]   = Sigma[[2]][2,1,4,2] = Sigma[[2]][1,2,2,4] = Sigma[[2]][2,1,2,4] = 1

        ##### Terceira derivada da matriz Sigma
        Sigma[[3]][1,1,4,2,2] = Sigma[[3]][1,1,2,4,2] = Sigma[[3]][1,1,2,2,4] = 2



        #####################################################################
        # Fim do caso lambda_e Correction = TRUE

      }else if(Correction == FALSE){
        # cat("Lambda_e - FALSE")

      }else{ stop("\n Declare 'Correction' sendo TRUE ou FALSE.") }

      ###################################################################
      # Comeco do caso lambda_x
    }else if(lambda_e == F & lambda_x > 0){
      # cat("Lambda_x")

      ##Estimativas de maxima verossimilhanca (e.m.v.)
      beta_hat     = ((lambda_x + 1)*S_xy)/(lambda_x*S_xx)
      mu_x_hat     = x_barra
      alfa_hat     = y_barra - (beta_hat*x_barra)
      sigma2_u_hat = S_xx/(lambda_x + 1)
      sigma2_e_hat = (lambda_x*(S_yy*S_xx - S_xy^2) - S_xy^2)/(lambda_x*S_xx)

      ##Estimativas de maxima verossimilhanca (e.m.v.) sob H0
      alfa_til     = y_barra - (beta_til *x_barra)
      mu_x_til     = x_barra;
      sigma2_u_til = S_xx/(lambda_x + 1);
      sigma2_e_til = ((lambda_x-1)*lambda_x*beta_til^2*S_xx - 2*(lambda_x+1)*lambda_x*beta_til*S_xy +	(lambda_x+1)^2*S_yy)/((lambda_x+1)^2)

      ##### Matrix de medias
      mu_hat[,] = matrix(c(alfa_hat + (beta_hat*mu_x_hat), mu_x_hat), ncol=1, nrow=2)
      mu_til[,] = matrix(c(alfa_til + (beta_til    *mu_x_til), mu_x_til), ncol=1, nrow=2)


      ##### Matrix de Covariancias
      Sigma_hat[,] = matrix(c((sigma2_e_hat)+(beta_hat^2*lambda_x*sigma2_u_hat),
                              beta_hat*lambda_x*(sigma2_u_hat),
                              beta_hat*lambda_x*(sigma2_u_hat),
                              ((lambda_x+1)*sigma2_u_hat)), ncol=2,nrow=2)


      Inv_Sigma_hat[,] = solve(Sigma_hat[,])
      Det_Sigma_hat= det(Sigma_hat[,])

      Sigma_til[,] = matrix(c((sigma2_e_til)+(beta_til^2*lambda_x*sigma2_u_til),
                              beta_til *lambda_x*(sigma2_u_til),
                              beta_til *lambda_x*(sigma2_u_til),
                              ((lambda_x+1)*sigma2_u_til)), ncol=2,nrow=2)

      Inv_Sigma_til[,] = solve(Sigma_til[,])
      Det_Sigma_til= det(Sigma_til[,])

      ##### Preenchendo os vetores Z, diferenca e produto
      for (i in 1:n)
      {
        Z[,,i]             = matrix(c(y[i],x[i]),ncol = 1,nrow = 2)

        diferenca_hat[,,i] = matrix(Z[,,i] - mu_hat[,],ncol = 1)
        diferenca_til[,,i] = matrix(Z[,,i] - mu_til[,],ncol = 1)

        Produto_hat[i]     = (t(matrix(diferenca_hat[,,i],ncol = 1)) %*% Inv_Sigma_hat[,] %*% matrix(diferenca_hat[,,i],ncol = 1))
        Produto_til[i]     = (t(matrix(diferenca_til[,,i],ncol = 1)) %*% Inv_Sigma_til[,] %*% matrix(diferenca_til[,,i],ncol = 1))
      }


      ###################################################################
      ###################################################################
      # Correction == TRUE - lambda_x
      if(Correction == TRUE){
        # cat("Lambda_x - TRUE")
        ###################################################################
        # Correction == TRUE - lambda_x

        ##### Primeira derivada do vetor mu
        Mu[[1]][1,1,1] = 1
        Mu[[1]][1,1,2] = mu_x_til
        Mu[[1]][1,1,3] = beta_til
        Mu[[1]][2,1,3] = 1


        ##### Segunda derivada do vetor mu
        Mu[[2]][1,1,2,3] = Mu[[2]][1,1,3,2] = 1

        ##### Primeira derivada da matriz Sigma
        Sigma[[1]][1,1,2] = 2*lambda_x*sigma2_u_til*beta_til
        Sigma[[1]][1,2,1] = Sigma[[1]][2,1,1] = lambda_x*sigma2_u_til
        Sigma[[1]][1,1,4] = 1
        Sigma[[1]][1,1,5] = beta_til^2*lambda_x
        Sigma[[1]][2,1,5] = Sigma[[1]][1,2,5] = beta_til*lambda_x
        Sigma[[1]][2,2,5] = lambda_x+1

        ##### Segunda derivada da matriz Sigma
        Sigma[[2]][1,1,2,2]   = 2*lambda_x*sigma2_u_til
        Sigma[[2]][1,1,5,2]   = Sigma[[2]][1,1,2,5] = 2*lambda_x*beta_til
        Sigma[[2]][1,2,5,2]   = Sigma[[2]][2,1,5,2] =   Sigma[[2]][1,2,2,5] = Sigma[[2]][2,1,2,5] = lambda_x

        ##### Terceira derivada da matriz Sigma
        Sigma[[3]][1,1,2,2,5] = Sigma[[3]][1,1,2,5,2] = Sigma[[3]][1,1,5,2,2] = 2*lambda_x

        #####################################################################
        # Fim do caso lambda_x Correction = TRUE


      }else if(Correction == FALSE){
        #cat("Lambda_x - FALSE")

      }else{ stop("\n Declare 'Correction' to be TRUE or FALSE. If TRUE, the output will print the corrected statistics and their p-values. \n") }
    }else{stop("\n Declare lambda_x or lambda_e ratios as a number greater than zero. \n Declare only one or the other. Use the one where you will admit to being known. \n")} # Fim da fucao if() de lambda_x e lambda_e


    ##### Log da verossimilhanca avaliados nas estimativas de maxima
    ##### verossimilhanca irrestrita e restrita

    L_teta_hat = - ((n/2)*log(Det_Sigma_hat)) - ((1/2)*sum(Produto_hat))
    L_teta_til = - ((n/2)*log(Det_Sigma_til)) - ((1/2)*sum(Produto_til))

    ##### Estatastica da razao de verossimilhancas
    LR = 2*(L_teta_hat - L_teta_til)

    ##### Valor de p
    p = 1
    VP_LR    = 1 -  stats::pchisq(LR, p )

    if(Correction){ # Segunda parte do correction TRUE
      ##### Primeira derivada da inversa de Sigma
      ##D_Inv_Sigma[[1]]
      for (r  in c(1:5))
      {D_Inv_Sigma[[1]][,,r] = -(Inv_Sigma_til[,] %*% Sigma[[1]][,,r] %*% Inv_Sigma_til[,])}


      ##### Segunda derivada da inversa de Sigma
      ##D_Inv_Sigma[[2]]
      for (r in 1:5)
      {
        for (s in 1:5)
        {
          D_Inv_Sigma[[2]][,,r,s] = (-(D_Inv_Sigma[[1]][,,r] %*% Sigma[[1]][,,s]   %*% Inv_Sigma_til[,])
                                     -(Inv_Sigma_til[,]      %*% Sigma[[2]][,,r,s] %*% Inv_Sigma_til[,])
                                     -(Inv_Sigma_til[,]      %*% Sigma[[1]][,,s]   %*% D_Inv_Sigma[[1]][,,r]))
        }
      }

      ##### Terceira derivada da inversa de Sigma
      ##D_Inv_Sigma[[3]]
      for (t in 1:5)
      {
        for (s in 1:5)
        {
          for (r in 1:5)
          {
            D_Inv_Sigma[[3]][,,r,s,t] = (-(D_Inv_Sigma[[2]][,,r,s]  %*% Sigma[[1]][,,t]    %*% Inv_Sigma_til[,])
                                         -(D_Inv_Sigma[[1]][,,s]    %*% Sigma[[2]][,,r,t]  %*% Inv_Sigma_til[,])
                                         -(D_Inv_Sigma[[1]][,,s]    %*% Sigma[[1]][,,t]    %*% D_Inv_Sigma[[1]][,,r])
                                         -(D_Inv_Sigma[[1]][,,r]    %*% Sigma[[2]][,,s,t]  %*% Inv_Sigma_til[,])
                                         -(Inv_Sigma_til[,]         %*% Sigma[[3]][,,r,s,t]%*% Inv_Sigma_til[,])
                                         -(Inv_Sigma_til[,]         %*% Sigma[[2]][,,s,t]  %*% D_Inv_Sigma[[1]][,,r])
                                         -(D_Inv_Sigma[[1]][,,r]    %*% Sigma[[1]][,,t]    %*% D_Inv_Sigma[[1]][,,s])
                                         -(Inv_Sigma_til[,]         %*% Sigma[[2]][,,r,t]  %*% D_Inv_Sigma[[1]][,,s])
                                         -(Inv_Sigma_til[,]         %*% Sigma[[1]][,,t]    %*% D_Inv_Sigma[[2]][,,r,s]))
          }
        }
      }


      ##### Quarta derivada da inversa de Sigma
      ##D_Inv_Sigma[[4]]
      for (r in 1:5)
      {
        for (s in 1:5)
        {
          for (t in 1:5)
          {
            for (w in 1:5)
            {
              D_Inv_Sigma[[4]][,,r,s,t,w]=(-(D_Inv_Sigma[[3]][,,r,s,t] %*% Sigma[[1]][,,w]       %*% Inv_Sigma_til[,])
                                           -(D_Inv_Sigma[[2]][,,s,t]   %*% Sigma[[2]][,,r,w]     %*% Inv_Sigma_til[,])
                                           -(D_Inv_Sigma[[2]][,,s,t]   %*% Sigma[[1]][,,w]       %*% D_Inv_Sigma[[1]][,,r])
                                           -(D_Inv_Sigma[[2]][,,r,t]   %*% Sigma[[2]][,,s,w]     %*% Inv_Sigma_til[,])
                                           -(D_Inv_Sigma[[1]][,,t]     %*% Sigma[[3]][,,r,s,w]   %*% Inv_Sigma_til[,])
                                           -(D_Inv_Sigma[[1]][,,t]     %*% Sigma[[2]][,,s,w]     %*% D_Inv_Sigma[[1]][,,r])
                                           -(D_Inv_Sigma[[2]][,,r,t]   %*% Sigma[[1]][,,w]       %*% D_Inv_Sigma[[1]][,,s])
                                           -(D_Inv_Sigma[[1]][,,t]     %*% Sigma[[2]][,,r,w]     %*% D_Inv_Sigma[[1]][,,s])
                                           -(D_Inv_Sigma[[1]][,,t]     %*% Sigma[[1]][,,w]       %*% D_Inv_Sigma[[2]][,,r,s])
                                           -(D_Inv_Sigma[[2]][,,r,s]   %*% Sigma[[2]][,,t,w]     %*% Inv_Sigma_til[,])
                                           -(D_Inv_Sigma[[1]][,,s]     %*% Sigma[[3]][,,r,t,w]   %*% Inv_Sigma_til[,])
                                           -(D_Inv_Sigma[[1]][,,s]     %*% Sigma[[2]][,,t,w]     %*% D_Inv_Sigma[[1]][,,r])
                                           -(D_Inv_Sigma[[1]][,,r]     %*% Sigma[[3]][,,s,t,w]   %*% Inv_Sigma_til[,])
                                           -(Inv_Sigma_til[,]          %*% Sigma[[3]][,,s,t,w]   %*% D_Inv_Sigma[[1]][,,r])
                                           -(D_Inv_Sigma[[1]][,,r]     %*% Sigma[[2]][,,t,w]     %*% D_Inv_Sigma[[1]][,,s])
                                           -(Inv_Sigma_til[,]          %*% Sigma[[3]][,,r,t,w]   %*% D_Inv_Sigma[[1]][,,s])
                                           -(Inv_Sigma_til[,]          %*% Sigma[[2]][,,t,w]     %*% D_Inv_Sigma[[2]][,,r,s])
                                           -(D_Inv_Sigma[[2]][,,r,s]   %*% Sigma[[1]][,,w]       %*% D_Inv_Sigma[[1]][,,t])
                                           -(D_Inv_Sigma[[1]][,,s]     %*% Sigma[[2]][,,r,w]     %*% D_Inv_Sigma[[1]][,,t])
                                           -(D_Inv_Sigma[[1]][,,s]     %*% Sigma[[1]][,,w]       %*% D_Inv_Sigma[[2]][,,r,t])
                                           -(D_Inv_Sigma[[1]][,,r]     %*% Sigma[[2]][,,s,w]     %*% D_Inv_Sigma[[1]][,,t])
                                           -(Inv_Sigma_til[,]          %*% Sigma[[3]][,,r,s,w]   %*% D_Inv_Sigma[[1]][,,t])
                                           -(Inv_Sigma_til[,]          %*% Sigma[[2]][,,s,w]     %*% D_Inv_Sigma[[2]][,,r,t])
                                           -(D_Inv_Sigma[[1]][,,r]     %*% Sigma[[1]][,,w]       %*% D_Inv_Sigma[[2]][,,s,t])
                                           -(Inv_Sigma_til[,]          %*% Sigma[[2]][,,r,w]     %*% D_Inv_Sigma[[2]][,,s,t])
                                           -(Inv_Sigma_til[,]          %*% Sigma[[1]][,,w]       %*% D_Inv_Sigma[[3]][,,r,s,t]))
            }
          }
        }
      }

      ### Funcao para calcular o traco de uma matriz
      traco = function(M){sum(diag(M))}

      for (t in 1:5)
      {
        for (w in 1:5)
        {
          ##### Matriz de Informacao de Fisher K

          K[t,w] = (n/2)*traco(D_Inv_Sigma[[1]][,,w] %*% Sigma[[1]][,,t]) - n*t(Mu[[1]][,,t]) %*% Inv_Sigma_til[,] %*% Mu[[1]][,,w]


          for (s in 1:5)
          {

            ### Q^s
            Q[t,s,w] = (
              -(n/2)*traco(
                D_Inv_Sigma[[2]][,,t,w] %*% Sigma[[1]][,,s] + D_Inv_Sigma[[2]][,,w,s]  %*% Sigma[[1]][,,t]
                + D_Inv_Sigma[[1]][,,w] %*% Sigma[[2]][,,t,s] + D_Inv_Sigma[[1]][,,t]  %*% Sigma[[2]][,,w,s]
                + Inv_Sigma_til[,]      %*% Sigma[[3]][,,t,w,s] + D_Inv_Sigma[[3]][,,t,w,s] %*% Sigma_til[,]
              )
              -n*
                (
                  + t(Mu[[1]][,,w])   %*% D_Inv_Sigma[[1]][,,t] %*% Mu[[1]][,,s]
                  + t(Mu[[1]][,,s])   %*% Inv_Sigma_til[,]      %*% Mu[[2]][,,t,w]
                  + t(Mu[[1]][,,w])   %*% Inv_Sigma_til[,]      %*% Mu[[2]][,,t,s]
                )
            )


            ### P^s
            P[t,s,w] = (
              -(n/2)*traco(
                D_Inv_Sigma[[2]][,,t,w]   %*% Sigma[[1]][,,s]
                + D_Inv_Sigma[[1]][,,w]     %*% Sigma[[2]][,,t,s]
                + D_Inv_Sigma[[1]][,,t]     %*% Sigma[[2]][,,w,s]
                + Inv_Sigma_til[,]          %*% Sigma[[3]][,,t,w,s]
                + D_Inv_Sigma[[3]][,,t,w,s] %*% Sigma_til[,]
              )
              -n*(
                t(Mu[[1]][,,t]) %*% D_Inv_Sigma[[1]][,,s] %*% Mu[[1]][,,w]
                + t(Mu[[1]][,,w]) %*% D_Inv_Sigma[[1]][,,t] %*% Mu[[1]][,,s]
                + t(Mu[[1]][,,t]) %*% D_Inv_Sigma[[1]][,,w] %*% Mu[[1]][,,s]
                + t(Mu[[1]][,,t]) %*% Inv_Sigma_til[,]      %*% Mu[[2]][,,w,s]
                + t(Mu[[1]][,,w]) %*% Inv_Sigma_til[,]      %*% Mu[[2]][,,t,s]
                + t(Mu[[1]][,,s]) %*% Inv_Sigma_til[,]      %*% Mu[[2]][,,t,w]
              )
            )


            for (r in 1:5)
            {
              ### A^{rs}
              A[r,w,t,s] = (
                -(n/8)*traco(
                  D_Inv_Sigma[[3]][,,t,w,r]   %*% Sigma[[1]][,,s]
                  + 4*D_Inv_Sigma[[3]][,,s,t,r]   %*% Sigma[[1]][,,w]
                  +   D_Inv_Sigma[[1]][,,t]       %*% Sigma[[3]][,,w,r,s]
                  +   D_Inv_Sigma[[1]][,,w]       %*% Sigma[[3]][,,t,r,s]
                  +   D_Inv_Sigma[[1]][,,r]       %*% Sigma[[3]][,,t,w,s]
                  +   D_Inv_Sigma[[2]][,,t,w]     %*% Sigma[[2]][,,r,s]
                  +   D_Inv_Sigma[[2]][,,w,r]     %*% Sigma[[2]][,,t,s]
                  + 5*D_Inv_Sigma[[2]][,,t,r]     %*% Sigma[[2]][,,w,s]
                  +   D_Inv_Sigma[[4]][,,t,w,r,s] %*% Sigma_til[,]
                )
                -(n/4)*(
                  t(Mu[[1]][,,t]) %*% D_Inv_Sigma[[2]][,,w,s] %*% Mu[[1]][,,r]
                  +   t(Mu[[1]][,,t]) %*% D_Inv_Sigma[[2]][,,w,r] %*% Mu[[1]][,,s]
                  - 3*t(Mu[[1]][,,t]) %*% D_Inv_Sigma[[2]][,,r,s] %*% Mu[[1]][,,w]
                  - 3*t(Mu[[1]][,,w]) %*% D_Inv_Sigma[[2]][,,t,s] %*% Mu[[1]][,,r]
                  +   t(Mu[[1]][,,s]) %*% D_Inv_Sigma[[2]][,,t,w] %*% Mu[[1]][,,r]
                  +   t(Mu[[1]][,,s]) %*% D_Inv_Sigma[[2]][,,t,r] %*% Mu[[1]][,,w]
                  +   t(Mu[[1]][,,t])	%*% D_Inv_Sigma[[1]][,,s]   %*%	Mu[[2]][,,w,r]
                  - 3*t(Mu[[1]][,,t])	%*% D_Inv_Sigma[[1]][,,r]   %*%	Mu[[2]][,,w,s]
                  +   t(Mu[[1]][,,t])	%*% D_Inv_Sigma[[1]][,,w]   %*%	Mu[[2]][,,r,s]
                  +   t(Mu[[1]][,,s])	%*% D_Inv_Sigma[[1]][,,t]   %*%	Mu[[2]][,,w,r]
                  +   t(Mu[[1]][,,s])	%*% D_Inv_Sigma[[1]][,,w]   %*%	Mu[[2]][,,t,r]
                  +   t(Mu[[1]][,,s])	%*% D_Inv_Sigma[[1]][,,r]   %*%	Mu[[2]][,,t,w]
                  +   t(Mu[[1]][,,r])	%*% D_Inv_Sigma[[1]][,,w]   %*%	Mu[[2]][,,t,s]
                  - 3*t(Mu[[1]][,,r])	%*% D_Inv_Sigma[[1]][,,t]   %*%	Mu[[2]][,,w,s]
                  +   t(Mu[[1]][,,r])	%*% D_Inv_Sigma[[1]][,,s]   %*%	Mu[[2]][,,t,w]
                  - 3*t(Mu[[1]][,,w])	%*% D_Inv_Sigma[[1]][,,s]   %*%	Mu[[2]][,,t,r]
                  - 3*t(Mu[[1]][,,w])	%*% D_Inv_Sigma[[1]][,,r]   %*%	Mu[[2]][,,t,s]
                  - 3*t(Mu[[1]][,,w])	%*% D_Inv_Sigma[[1]][,,t]   %*%	Mu[[2]][,,r,s]
                  +   t(Mu[[2]][,,t,w]) %*% Inv_Sigma_til[,]      %*% Mu[[2]][,,r,s]
                  - 3*t(Mu[[2]][,,t,r]) %*% Inv_Sigma_til[,]      %*% Mu[[2]][,,w,s]
                  +   t(Mu[[2]][,,w,r]) %*% Inv_Sigma_til[,]      %*% Mu[[2]][,,t,s]
                )
              )

            }## Fim do loop do r
          }## Fim do loop do s
        }## Fim do loop do w
      } ## Fim do loop do t

      Inv_K[,] = (solve(K[,]))


      for (r in 1:5)
      {
        for (s in 1:5)
        {
          ##### L
          L[r,s] = traco(Inv_K[,] %*% A[,,r,s])

          ##### M
          M[r,s] = (
            -(1/6)*traco(Inv_K[,] %*% P[,,r] %*% Inv_K[,] %*%   P[,,s])
            +traco(Inv_K[,] %*% P[,,r] %*% Inv_K[,] %*% t(Q[,,s]))
            -traco(Inv_K[,] %*% Q[,,r] %*% Inv_K[,] %*%   Q[,,s])
          )

          ##### N
          N[r,s] = (
            -(1/4)*traco(P[,,r] %*% Inv_K[,])*traco(P[,,s] %*% Inv_K[,])
            +traco(P[,,r] %*% Inv_K[,])*traco(Q[,,s] %*% Inv_K[,])
            -traco(Q[,,r] %*% Inv_K[,])*traco(Q[,,s] %*% Inv_K[,])
          )
        }
      }


      ##### Epsilon(p)
      Epsilon_p = traco(Inv_K[,] %*% (L[,] - M[,] - N[,]))

      ##### Calculo de Epsilon(p-1):
      ##### Retiramos a segunda coluna e a segunda linha das matrizes, pois teta = (alfa, beta,vmu_x, sigma2_x, sigma2_u)
      a = 2
      K_p_1[,] = K[-a,-a]
      L_p_1[,] = L[-a,-a]
      M_p_1[,] = M[-a,-a]
      N_p_1[,] = N[-a,-a]

      Inv_K_p_1[,] = solve(K_p_1[,])

      Epsilon_p_1 = traco(Inv_K_p_1[,] %*% (L_p_1[,] - M_p_1[,] - N_p_1[,]))

      ##### Constante C (Correcao de Bartlett)
      CB1 = 1 + (Epsilon_p - Epsilon_p_1)
      CB2 = exp(1 - CB1)
      CB3 = 1 - (Epsilon_p - Epsilon_p_1)

      ##### Estatastica da razao de verossimilhancas corrigidas por Bartlett
      LR_B1 = LR/CB1
      LR_B2 = LR*CB2
      LR_B3 = LR*CB3

      ##### Valor de p
      VP_LR_B1 = 1 - stats::pchisq(LR_B1, p )
      VP_LR_B2 = 1 - stats::pchisq(LR_B2, p )
      VP_LR_B3 = 1 - stats::pchisq(LR_B3, p )

    } # Fim da segunda parte da correcao

    ##################################################################
    ##### Impressao de resultados




    # Function do texto de saida
    if(Correction){
      ##### Impressão de resultados - Correction TRUE
      if(lambda_e == F & lambda_x > 0){ # Caso lamda_x

        Coefficients = data.frame( "."  = c("(Intercept)","x"),
                                   Estimate      = round( c(alfa_hat, beta_hat),4),
                                   "Estimate_H0" = round( c(alfa_til, beta_til),4),
                                   LR            = c(round( LR,4),"---"),
                                   "P_Value_LR"  = c(format(VP_LR, scientific = T, digits = 3),if(VP_LR < 0.001){"***"}else if(VP_LR < 0.01){"**"}else if(VP_LR < 0.05){"*"}else if(VP_LR < 0.1){"."}else{""}) )

        dados = data.frame("Statistics"= c("LR", 'LR_B1', 'LR_B2', 'LR_B3'),
                           "Value"   = round( c(LR, LR_B1, LR_B2, LR_B3),4),
                           "P_Value" = c(format(VP_LR, scientific = T, digits = 3),
                                         format(VP_LR_B1, scientific = T, digits = 3),
                                         format(VP_LR_B2, scientific = T, digits = 3),
                                         format(VP_LR_B3, scientific = T, digits = 3)  ),
                           "."        = c(if(VP_LR < 0.001){"***"}else if(VP_LR < 0.01){"**"}else if(VP_LR < 0.05){"*"}else if(VP_LR < 0.1){"."}else{""},
                                          if(VP_LR_B1 < 0.001){"***"}else if(VP_LR_B1 < 0.01){"**"}else if(VP_LR_B1 < 0.05){"*"}else if(VP_LR_B1 < 0.1){"."}else{""},
                                          if(VP_LR_B1 < 0.001){"***"}else if(VP_LR_B2 < 0.01){"**"}else if(VP_LR_B2 < 0.05){"*"}else if(VP_LR_B2 < 0.1){"."}else{""},
                                          if(VP_LR_B1 < 0.001){"***"}else if(VP_LR_B3 < 0.01){"**"}else if(VP_LR_B3 < 0.05){"*"}else if(VP_LR_B3 < 0.1){"."}else{""} ))



        # Inicio da funcao do texto de saida
        Saida = function(){
          cat("\n")
          cat("Call: \n")
          cat("MEM( formula = y ~ x, = lambda_x = ", lambda_x, " Correction = TRUE)" ," \n")
          cat("\n")
          cat("Coefficients: \n")
          print.data.frame(Coefficients,row.names = F)
          cat("\n")
          cat("Bartlett's correction: \n")
          print.data.frame(dados,row.names = F)
          cat("--- \n")
          cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 \n")
          cat("\n")
          # R2 pagina 7
          r2 = (S_xy^2)/(S_xx*S_yy)
          cat("Multiple R-squared:", round(r2,3), ", Mean of x:", x_barra, ", Mean of y:", y_barra, "\n")
          cat("Sigma_e: ", round(sigma2_e_hat,3), "considering the null hypothesis:",  round(sigma2_e_til,3),"\n")
          cat( "Sigma_u: ", round(sigma2_u_hat,3), "considering the null hypothesis:",  round(sigma2_u_til,3))
          cat("\n \n")
        } #End of function exit
      }else if(lambda_x == F & lambda_e > 0){ #Correction TRUE lambd_e
        Coefficients = data.frame( "."  = c("(Intercept)","x"),
                                   Estimate      = round( c(alfa_hat, beta_hat),4),
                                   "Estimate_H0" = round( c(alfa_til, beta_til),4),
                                   LR            = c(round( LR,4),"---"),
                                   "P_Value_LR"  = c(format(VP_LR, scientific = T, digits = 3),if(VP_LR < 0.001){"***"}else if(VP_LR < 0.01){"**"}else if(VP_LR < 0.05){"*"}else if(VP_LR < 0.1){"."}else{""}) )

        dados = data.frame("Statistics"= c("LR", 'LR_B1', 'LR_B2', 'LR_B3'),
                           "Value"   = round( c(LR, LR_B1, LR_B2, LR_B3),4),
                           "P_Value" = c(format(VP_LR, scientific = T, digits = 3),
                                         format(VP_LR_B1, scientific = T, digits = 3),
                                         format(VP_LR_B2, scientific = T, digits = 3),
                                         format(VP_LR_B3, scientific = T, digits = 3)  ),
                           "."        = c(if(VP_LR < 0.001){"***"}else if(VP_LR < 0.01){"**"}else if(VP_LR < 0.05){"*"}else if(VP_LR < 0.1){"."}else{""},
                                          if(VP_LR_B1 < 0.001){"***"}else if(VP_LR_B1 < 0.01){"**"}else if(VP_LR_B1 < 0.05){"*"}else if(VP_LR_B1 < 0.1){"."}else{""},
                                          if(VP_LR_B1 < 0.001){"***"}else if(VP_LR_B2 < 0.01){"**"}else if(VP_LR_B2 < 0.05){"*"}else if(VP_LR_B2 < 0.1){"."}else{""},
                                          if(VP_LR_B1 < 0.001){"***"}else if(VP_LR_B3 < 0.01){"**"}else if(VP_LR_B3 < 0.05){"*"}else if(VP_LR_B3 < 0.1){"."}else{""} ))



        # Inicio da funcao do texto de saida
        Saida = function(){
          cat("\n")
          cat("Call: \n")
          cat("MEM( formula = y ~ x, = lambda_e = ", lambda_e, " Correction = TRUE)" ," \n")
          cat("\n")
          cat("Coefficients: \n")
          print.data.frame(Coefficients,row.names = F)
          cat("\n")
          cat("Bartlett's correction: \n")
          print.data.frame(dados,row.names = F)
          cat("--- \n")
          cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 \n")
          cat("\n")
          # R2 pagina 7
          r2 = (S_xy^2)/(S_xx*S_yy)
          cat("Multiple R-squared:", round(r2,3), ", Mean of x:", x_barra, ", Mean of y:", y_barra, "\n")
          cat("Sigma_x: ", round(sigma2_x_hat,3), "considering the null hypothesis:",  round(sigma2_x_til,3),"\n")
          cat( "Sigma_u: ", round(sigma2_u_hat,3), "considering the null hypothesis:",  round(sigma2_u_til,3))
          cat("\n \n")
        } #End of function exit
      } # FIm do caso Correction TRUE lambd_e
    }else{
      ##### Impressão de resultados - Correction FALSE
      if(lambda_e == F & lambda_x > 0){ # Caso lamda_x
        Coefficients = data.frame( "."  = c("(Intercept)","x"),
                                   Estimate      = round( c(alfa_hat, beta_hat),4),
                                   "Estimate_H0" = round( c(alfa_til, beta_til),4),
                                   LR            = c(round( LR,4),"---"),
                                   "P_Value_LR"  = c(format(VP_LR, scientific = T, digits = 3),if(VP_LR < 0.001){"***"}else if(VP_LR < 0.01){"**"}else if(VP_LR < 0.05){"*"}else if(VP_LR < 0.1){"."}else{""}) )

        # Inicio da funcao do texto de saida
        Saida = function(){
          cat("\n")
          cat("Call: \n")
          cat("MEM( formula = y ~ x, = lambda_x = ", lambda_x, " Correction = TRUE)" ," \n")
          cat("\n")
          cat("Coefficients: \n")
          print.data.frame(Coefficients,row.names = F)
          cat("--- \n")
          cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 \n")
          cat("\n")
          # R2 pagina 7
          r2 = (S_xy^2)/(S_xx*S_yy)
          cat("Multiple R-squared:", round(r2,3), ", Mean of x:", x_barra, ", Mean of y:", y_barra, "\n")
          cat("Sigma_e: ", round(sigma2_e_hat,3), "considering the null hypothesis:",  round(sigma2_e_til,3),"\n")
          cat( "Sigma_u: ", round(sigma2_u_hat,3), "considering the null hypothesis:",  round(sigma2_u_til,3))
          cat("\n \n")
        } #End of function exit
      }else if(lambda_x == F & lambda_e > 0){ #Correction FALSE lambd_e
        Coefficients = data.frame( "."  = c("(Intercept)","x"),
                                   Estimate      = round( c(alfa_hat, beta_hat),4),
                                   "Estimate_H0" = round( c(alfa_til, beta_til),4),
                                   LR            = c(round( LR,4),"---"),
                                   "P_Value_LR"  = c(format(VP_LR, scientific = T, digits = 3),if(VP_LR < 0.001){"***"}else if(VP_LR < 0.01){"**"}else if(VP_LR < 0.05){"*"}else if(VP_LR < 0.1){"."}else{""}) )

        # Inicio da funcao do texto de saida
        Saida = function(){
          cat("\n")
          cat("Call: \n")
          cat("MEM( formula = y ~ x, = lambda_e = ", lambda_e, " Correction = TRUE)" ," \n")
          cat("\n")
          cat("Coefficients: \n")
          print.data.frame(Coefficients,row.names = F)
          cat("--- \n")
          cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 \n")
          cat("\n")
          # R2 pagina 7
          r2 = (S_xy^2)/(S_xx*S_yy)
          cat("Multiple R-squared:", round(r2,3), ", Mean of x:", x_barra, ", Mean of y:", y_barra, "\n")
          cat("Sigma_x: ", round(sigma2_x_hat,3), "considering the null hypothesis:",  round(sigma2_x_til,3),"\n")
          cat( "Sigma_u: ", round(sigma2_u_hat,3), "considering the null hypothesis:",  round(sigma2_u_til,3))
          cat("\n \n")
        } #End of function exit

      }
    } # Fim do if de print da saida em Correction TRUE ou FALSE

    ##########################################################################
    # Saida com os resultados
    return(Saida())

  }else{ stop("\n INVALID PARAMETER IN MEM() FUNCTION \n The following conditions must be respected: \n - Declare x and y as a vector of numbers of the same size. \n - Declare lambda_x or lambda_e} reasons as a number greater than zero, declare only one or the other. Use the one where you will admit to being known. \n - Declare beta_til as a number, remembering that H_0: beta = beta_til. \n - Declare 'Correction' to be TRUE or FALSE. If TRUE, the output will print the corrected statistics and their p-values. \n", noBreaks. = TRUE) } # Fim da funcao geral If()
} # Fim da funcao MEM
