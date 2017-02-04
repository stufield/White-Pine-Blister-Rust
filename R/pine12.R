##################################################
### 12 CLASS INFECTION MODEL #####################
### without Genetics 12x12 #######################
##################################################

#' The pine12 Model
#' 
#' The full 6 class with susceptible & infected in the pine12 model. There are
#' three forms of non-linearity: 2 forms of density-dependent seed germination
#' (rALs & r.cache) & density-dependent fecundity (r.cones). Genetics or
#' resistance is added to the model in \code{pine36}.
#' 
#' State.vector = c(SD, SD1, SD2, SP, YA, MA, SDi, SD1i, SD2i, SPi, YAi, MAi)
#' Beta for the Adult stage (Range = 0.016 - 0.14); default = 0.044.
#' 
#' @param Gen Integer. Number of generations to project the population.
#' @param x1 A vector indicating the initial population.
#' @param Beta The transmission probability
#' @param LAIb Background LAI due to external species in the simulation.
#' @param r.site Suitability of the site to germination. Defaults to no effect (r.site=1).
#' @param delta Effect of infection on survivorship.
#' @param C.f Effect of infection on fecundity.
#' @param rho Relative difference in fecundity between young adults (jv) &
#' mature adults (ad).
#' @param Cmax Maximum cone production per reproductive tree.
#' @param SpC Seeds per cache
#' @param S.cone Seeds per cone
#' @param P.find Proportion of seeds found by dispersal vectors (birds).
#' @param nBirds Number of birds per plot.
#' @param P.cons Proportion of seeds immediately consumed by birds.
#' @param alpha1 Coefficient corresponding to the LAI of the sapling (SA)
#' stage.
#' @param alpha2 Coefficient for the conversion of dbh -> LAI.
#' @param alpha3 Coefficient in the exponent for the conversion of dbh -> LAI.
#' @param c2 Cost of infection to survivorship of SD1 individuals.
#' @param c3 Cost of infection to survivorship of SD2 individuals.
#' @param dbh.v A vector (of length=6) describing the mean dbh for each class
#' (1 -> 6).
#' @param M Numeric. A vector (of length=5) describing the mortality for each class (1
#' -> 5). Subsequently converted to survivorship by the \code{vitals}
#' functions.
#' @param m6 The mortality of the adult class (6).
#' @param R Numeric. A vector (of length=6) describing the mean Residence time for each
#' class (1 -> 6).
#' @param plot Logical. Should output plots be created?
#' @param csv Logical. Should all solutions and time steps be converted to a
#' 2-dim matrix and exported to a *.csv file (in the working directory)?
#' @param filename Character string. If \code{csv=TRUE}, the name of the csv
#' file to be produced.
#' @param silent Logical. Should the output be silenced?
#' @return A list containing:
#' \item{theta }{List of model parameters}
#' \item{InitialPop }{A vector of the initial population (x1 -> x12)}
#' \item{PopProjection }{A matrix of the population projection}
#' \item{PopProjGroupClass }{A matrix of the population projection grouped 
#' by class}
#' \item{PopTotals }{A vector of the total population by over each 
#' time step in the projection}
#' \item{Prevalence }{A matrix of the infection prevalence (proportion infected) 
#' at each time step in the projection, grouped by class & of the total population}
#' \item{LAI.y }{A vector of LAI of the intermediate population, following
#' survivorship & transition but prior to fecundity & infection}
#' \item{r.cache }{A vector of the germination coefficient at each time step in the
#' projection} 
#' \item{ConesTree }{A vector of the cones per tree (density-dependent) at 
#' each time step in the projection}
#' \item{GerminationRate }{A vector of the germination rate
#' (density-dependent) at each time step in the projection}
#' \item{Smatrix }{A matrix of the basic projection matrix (minus fecundity)}
#' \item{Bmatrix }{A matrix of the infection matrix}
#' \item{y.vec }{A vector of the final intermediate population in the projection
#' corresponding to y.vec_(Gen)}
#' \item{f.vec }{A vector of the final fecundity with non-zero entries in the
#' entries 1 & 2 only. Entry 1 is the new seeds, entry 2 is the new SD1
#' individuals, corresponding to seed production and germination respectively}
#' \item{FinalPopVec }{A vector of the final solution of the projection}
#' \item{FinalSum }{The sum of the final solution of the projection}
#' \item{Adults }{The sum of the adults (MA) only in the final time step}
#' \item{DeadTrees }{A matrix of the cumulate sum of dead trees by class and
#' time step}
#' \item{LambdaGrowth }{A vector of the proportional growth in successional time
#' steps during the projection}
#' @author Stu Field
#' @seealso \code{\link{pine36}}
#' @keywords pine12 pine36 CSU white pine blister rust
#' @examples
#' 
#' pine12(Gen=25, silent=FALSE)
#' 
#' # Eqm solution for LAIb = 2.0
#' ipop <- c(35504, 22, 45, 36, 51, 201, 0, 0, 0, 0, 0, 0)
#' pine12(x1=ipop)
#'
#' # Eqm solution
#' \dontrun{
#' pine12(Gen=2500)$FinalPopVec
#' }
#' 
#' @importFrom stuRpkg diagR zeros BlockMat
#' @importFrom utils write.csv
#' @importFrom graphics grid legend lines par
#' @export pine12
pine12 <- function(Gen = 25,
                   x1 = c(62580, 38, 79, 65, 91, 353, 0,0,0,0,0,0),
                   #x1 = c(35504, 22, 45, 36, 51, 201, 0,0,0,0,0,0), # LAIb=2 eqm soln
                   Beta = 0.044,
                   LAIb = 0,
                   r.site = 1,
                   delta = 0.15,
                   C.f = 1/8,
                   rho =0.1,
                   Cmax = 7.5,
                   SpC = 3.7,
                   S.cone = 46,
                   P.find = 0.8,
                   nBirds = 3,
                   P.cons = 0.3,
                   alpha1 = 0.456, 
                   alpha2 = 0.0736,
                   alpha3 = 2.070,
                   c2 = 0.99, 
                   c3 = 0.87,
                   dbh.v = c(0, 0, 0, 2.05, 12.5, 37),
                   M = c(1, 0.152, 0.105, 0.02, 0.015),
                   m6 = 0.005,
                   R = c(1, 4, 16, 20, 50, 0),
                   plot = FALSE,
                   csv = FALSE,
                   filename = NULL,
                   silent = TRUE) {      

   #####################################################
   # Create storage matrices & vectors
   #####################################################
   n_stages <- length(x1)

   # overall storage of population projection
   popn.Mat <- matrix(0, Gen, n_stages,
                      dimnames=list(paste0("Gen",1:Gen), paste0("Stage", 1:n_stages)))

   popn.Mat[1,] = x1                             # first row is the initial population
   Pop.Total    = sum(x1[-1])                    # Remove seeds from total population
   # Initiate storage vectors
   LAI.v        = numeric(Gen)                   # LAI
   Prev.v       = sum(x1[8:12]) / (sum(x1[-1]))  # prevalence
   r.cache.v    = numeric(Gen)                   # r.cache
   Cones.v      = numeric(Gen)                   # Cones per tree
   SR.v         = numeric(Gen)                   # Seedling Recruitment

   ########################################
   #### Calculate Vital Rates #############
   ########################################
   M <- c(M, m6)
   S.par <- vitals(R, M)["surv", ]
   T.par <- vitals(R, M)["trans", ]

   ########################################
   #### Set linear Parameters #############
   ########################################
   ### Leaf Area Calculation ##############
   ########################################
   LA.v    <- alpha2 * (dbh.v^alpha3)
   LA.v[3] <- alpha1     # Don't forget the sedondary seedlings

   ##########################
   #### Beta parameters #####
   ##########################
   B.par <- c(0, rep(Beta, 5))       # Beta is constant for all classes; Used to make Beta Matrix (BM)
   #B.par <- (LA.v / LA.v) * Beta    # This is if Beta is class-dependent

   ##########################
   #### Cost parameters #####
   ##########################
   C.par <- 1 - exp(-delta * (dbh.v))   # Cost of infection to survivorship
   C.par[2] <- 1 - c2
   C.par[3] <- 1 - c3
   C.par <- c(rep(1, 6), C.par)         # Make into 12 entry matrix, 1s for uninfected classes

   ####################################
   #### Set up projection matrix ######
   ######## Linear Map Matrix #########
   ####################################
   # Survivorship Matrix
   S <- diag(S.par) + diagR(T.par[1:5], -1)

   # Cost Matrix
   C <- diag(C.par)

   # Beta Matrix
   B <- diag(B.par)

   # 6x6 Identity
   I6 <- diag(6)

   ### Combine sub-matrices using MatComb function
   SM <- BlockMat(list(S, zeros(6), zeros(6), S), b=2) %*% C

   # Combine sub-matrices using MatComb function
   BM <- BlockMat(list(I6-B, zeros(6), B, I6), b=2)

   
   #  Loop over generations #
   for ( n in 2:Gen ) {
      ######################################
      # Reset for new iteration
      ######################################
      f1.vec <- rep(0, n_stages)
      f2.vec <- rep(0, n_stages)
      x.vec  <- popn.Mat[n-1,]
      ################################
      # LAI.x Calculation
      ################################
      # Projected LAI of last year's popn (x.vec)
      # 1 ha = 10000 m^2
      LAI.x <- sum(LA.v * x.vec) / 10000 + LAIb   # additiion of background LAI

      #####################
      # Determine f[2]
      #####################
      SpB <- x.vec[1] / nBirds
      r.cache <- 0.73 / (1 + exp((31000 - SpB)/3000)) + 0.27
      r.cache.v[n] <- r.cache               # Store & track r.cache
      r.ALs <- 1 / (1 + exp(2 * (LAI.x - 3)))
   
      ############################
      # Seedling Recruitment
      # proportion seeds becoming seedlings
      ############################
      SR <- (((1 - P.find) * (1 - P.cons)) / SpC) * r.cache * r.ALs * r.site   # r.site
      SR.v[n] <- SR                    # Store & track SR for plotting
      f2.vec[2] <- x.vec[1] * SR       # Determine number of seedlings f_2
      x.vec <- x.vec + f2.vec          # Add seedlings to population
   
      ############################
      # Survive & transition
      ############################
      y.vec <- SM %*% x.vec            # Matrix multiplication of population trajectory

      ################################
      # LAI.y Calculation
      ################################
      LAI.y <- sum(LA.v * y.vec) / 10000 + LAIb   # LAI of this year's popn (sgf 12/16)
      LAI.v[n] <- LAI.y                           # Store and track of LAI over time

      #####################
      # Determine f[1]
      #####################
      r.cones    <- (0.5/(1 + exp(5*(LAI.y - 2.25))) + 0.5)
      C.tree     <- r.cones * Cmax
      Cones.v[n] <- C.tree
      f1.vec[1]  <- (S.cone * C.tree) * (rho*y.vec[5] + y.vec[6] + C.f * (rho*y.vec[11] + y.vec[12]))
   
      ######################################################################
      # add f_1(y.vec) the nonlinear function to y.vec & infect with B
      ######################################################################
      y.vec <- y.vec + f1.vec         # Add seeds to population; no seed bank, thus S[1,1]=0 & y.vec is 0 until here
      popn.Mat[n,] <- BM %*% y.vec    # Infection process to population in Fall

      ###############################################
      # Calculate some variables to follow
      ###############################################
      Pop.Total[n] <- sum(popn.Mat[n, -1])                             # Remove seedlings from population count
      Prev.v[n]    <- sum(popn.Mat[n, 8:12]) / (sum(popn.Mat[n, -1]))   # Calculate prevalence of rust in popn excluding seeds

   }
   # END OF LOOP

   theta <- list(Gen = Gen,
                 Transmission = Beta,
                 LAIb = LAIb, 
                 r.site = r.site,
                 delta = delta, 
                 C.f = C.f, 
                 rho = rho,
                 Cmax = Cmax, 
                 SpC = SpC, 
                 S.cone = S.cone,
                 P.find = P.find, 
                 nBirds = nBirds, 
                 P.cons = P.cons,
                 alpha1 = alpha1, 
                 alpha2 = alpha2,
                 alpha3 = alpha3,
                 c2 = c2,
                 c3 = c3,
                 plot = plot,
                 silent = silent
            )

   theta$tree_pars <- data.frame(rbind(dbh.v=dbh.v, M=M, R=R))
   names(theta$tree_pars) <- paste0("Stage", seq_along(dbh.v))

   #############################
   # Post-process variables
   #############################
   Dead      <- calcDeadTrees(popn.Mat, s=S.par, t=T.par, cc=C.par)   # calculate dead trees per class
   Groups    <- groupClass(popn.Mat)      # group projection by susc & inf classes
   PrevGroup <- prevClass(popn.Mat, Prev.v)
   Lambda    <- LambdaGrow(Pop.Total)

   #########################################
   # Convert Solutions & Export to CSV
   #########################################
   if ( csv ) {
      if ( is.null(filename) )
         stop("Please pass a file name for the CSV output file via the filename= argument")
      rownames(popn.Mat) <- paste0("t", 1:Gen)
      colnames(popn.Mat) <- paste0("Class_", 1:n_stages)
      write.csv(popn.Mat, file=sprintf("%s_FullSoln_%s.csv", filename, Sys.Date()))
   }

   ############################
   # Plot variables of interest
   # one day pull this plotting routine out into separate function
   ############################
   if ( plot ) {
      par(mfrow=c(3, 3))
      plot(1:Gen, Pop.Total, type="l", xlim=c(0, Gen), xlab="Years",
           main="Total Population", 
           ylab="Total Whitebark Population: Class 2-12", 
           col="navy", lwd=1.5)
      grid()
      plot(1:Gen, popn.Mat[,1], type="l", xlim=c(0,Gen), xlab="Years", 
           main="Seed Population", ylab="Total Seeds", col="darkgreen",
           lwd=1.5, lty=1)
      grid()
      plot(1:Gen, PrevGroup[, "total"], type="l", col="darkred", xlab="Years",
           lwd=1.5, main="Rust Prevalence",
           ylab="Proportion infected individuals")
      grid()

      for ( p in 2:n_stages ) {
         if ( p==2 ) {
            plot(1:Gen, popn.Mat[,p], type="l",
                 ylim=c(0, max(popn.Mat[, 2:n_stages])),
                 xlim=c(0, Gen), xlab="Years",
                 main="Susceptible Population",
                 ylab="Whitebark Population by Class", col=p, lwd=1.5)
            grid()
            legend("topright", legend=c(2:6), lty=c(rep(1,5)), col=c(2:6), 
                   lwd=1.5, cex=0.4, bg="gray95")
         }

         if ( p >= 3 && p <= 6 )
            lines(1:Gen, popn.Mat[, p], type="l", xlim=c(0, Gen), col=p, lwd=1.5)

         if ( p == 8 ) {
            plot(1:Gen, popn.Mat[,p], type="l", 
                 ylim=c(0,max(popn.Mat[, 9:n_stages])), 
                 xlim=c(0,Gen),
                 xlab="Years", 
                 main="Infected Population", 
                 ylab="Whitebark Population by Class", 
                 lty=2, col=p-6, lwd=1.5)
            grid()
            legend("topright", legend=c(8:12), lty=c(rep(2,5)),
                   col=c(2:6), lwd=1.5, cex=0.4, bg="gray95")
         }
         if ( p >= 9 )
            lines(1:Gen, popn.Mat[,p], type="l", xlim=c(0,Gen), lty=2,
                  col=p-6, lwd=1.5)
      }

      plot(1:(Gen-1), LAI.v[2:Gen], type="l", xlab="Years",
           main="Leaf Area Index", ylab="LAI.y", lty=4, lwd=1.5)
      grid()
      plot(1:(Gen-1), SR.v[2:Gen], type="l", xlim=c(0,Gen), xlab="Years", 
           main="Seedling Recruitment", ylab="Germination rate", 
           col="darkorchid", lwd=1.5, lty=4)
      grid()
      plot(1:(Gen-1), r.cache.v[2:Gen], type="l", xlab="Years", 
              ylab="r.cache", main="r.cache", col="red", lwd=1.5, lty=4)
      grid()
      plot(1:(Gen-1), Cones.v[2:Gen], type="l", xlab="Years", ylab="C.tree", 
           main="Total Cones per Tree", sub="A function of r.cones",
           col="brown", lwd=1.5, lty=4)
      grid()
   }
      
   ret <- list(theta = theta,
               InitialPop = x1,
               PopProjection = round(popn.Mat, 3),
               PopProjGroupClass = round(Groups, 3),
               PopTotals = round(Pop.Total, 2),
               Prevalence = PrevGroup,
               LAI.y = LAI.v,
               r.cache = r.cache.v,
               Cones.tree = Cones.v,
               GerminationRate = SR.v,
               Smatrix = SM,
               Bmatrix = BM,
               y.vec = y.vec,
               f.vec = f1.vec + f2.vec,
               FinalPopVec = popn.Mat[Gen, ],
               FinalSum = sum(popn.Mat[Gen, -1]),
               Adults = sum(popn.Mat[Gen, c(6, 12)]),
               DeadTrees = Dead,
               LambdaGrowth = Lambda
         )

   if ( silent )
      invisible(ret)
   else
      ret
}


