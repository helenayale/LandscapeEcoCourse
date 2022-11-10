#functions for simulating eco-evolutionary dynamics of metapopulations
initialize<-function(n.patches,n.init)
{
  adults <- data.frame(patch=rep(1:n.patches,each=n.init)) #each patch receives n.init individuals
  adults$a1 <- sample(100,size=nrow(adults),replace=T)/100 #diploid genotypes of these individuals 
  adults$a2 <- sample(100,size=nrow(adults),replace=T)/100 #are drawn randomly
  return(adults) #return the initial adults data frame
}

dispersal<-function(adults,d,n.patches)
{
  #sample the individuals that will disperse
  dispersing <- sample(nrow(adults),size=nrow(adults)*d) 
  #for each of the dispersing individuals, sample the patch to which they disperse
  adults[dispersing,"patch"]<-sample(n.patches,size=length(dispersing),replace=T)
  #return the updated data.frame
  return(adults)
}

fitness <- function(adults,neutral,w.max,env=NULL)
{
  if (neutral) adults$w<-w.max #in the neutral case all indviduals have fitness=w.max
  else #non-neutral case
  {
    adults$b<-adults$a1+adults$a2 #calculate the breeding value of all individuals
    adults$w<- w.max*exp(-(env[adults$patch]-adults$b)^2) #the fitness of each individual depends on how close its breeding value matches the local environment
  }
  #return the updated data.frame
  return(adults)
}

reproduction <- function(adults,k)
{
  offspring<-adults[NULL,] #create an empty data frame that will contain the offspring
  adults$n.offspring <- rpois(n=nrow(adults),lambda=adults$w) #for each adult, the number of offspring is drawn from a Poisson distribution with mean equal to the individual's fitness
  adults<-adults[adults$n.offspring>0,] #delete all adults that do not have offspring
  if (nrow(adults)>0) #if there are still adults left..
  {
    for (p in unique(adults$patch)) #..loop through all populations..
    {
      pop<-adults[adults$patch==p,] #extract the reproducing adults in the local population
      gametes<-c(rep(pop$a1,times=pop$n.offspring),rep(pop$a2,times=pop$n.offspring)) #create the gamete pool
      n.recruits<-min(k,sum(pop$n.offspring)) #competition: the number of successfully recruiting offspring cannot exceed the carrying capacity K
      recruit.gametes<-sample(gametes,size=2*n.recruits,replace=F) #random mating: the diploid genotypes of the recruits are sampled randomly from the gamete pool
      offspring<-rbind(offspring,data.frame(patch=p,a1=recruit.gametes[1:n.recruits],
                                            a2=recruit.gametes[(n.recruits+1):(2*n.recruits)])) #assemble the data for the new recruits in the population and add it to the offspring data frame
    } 
  } else offspring <- adults #if there are no adults left, the offspring are the adults (that is there are no offspring either)
  return(offspring) #return the offspring data frame
}



