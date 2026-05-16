

## Input: d
d=2
n=10

u=runif(d, -1,1)
sum(u^2)
rm(u)

set.seed(1)
resultados = matrix(NA, nrow = n, ncol = d)
for(i in 1:n){
  norma2=2
  while(norma2>1){
    u=runif(d, -1,1)
    if(sum(u^2)<=1) norma2 = sum(u^2) 
  }
  resultados[i,] <- u/sqrt(norma2)
}
resultados


rspherical = function(d, n){
  resultados = matrix(NA, nrow = n, ncol = d)
  for(i in 1:n){
    norma2=2
    while(norma2>1){
      u=runif(d, -1,1)
      if(sum(u^2)<=1) norma2 = sum(u^2) 
    }
    resultados[i,] <- u/sqrt(norma2)
  }
  resultados
}



x11()
set.seed(1)
plot(rspherical(d=2, n=100))


library(plotly)

# Criar uma matriz para armazenar densidades
z <- rspherical(d=3, n=100)

?plot_ly

# Criar gráfico interativo com plotly
plot_ly(data = data.frame(z),
  # x = z[,1],
  # y = z[,2],
  # z = z[,3],
  type = "scatter",
  colorscale = "Viridis"
) 
#%>%
  # layout(
  #   title = list(text = "AAAAA", font = list(size = 16)),
  #   scene = list(
  #     xaxis = list(title = "x"),
  #     yaxis = list(title = "y"),
  #     zaxis = list(title = "z")
  #   )
  # )



# Input: mu, Sigma



























