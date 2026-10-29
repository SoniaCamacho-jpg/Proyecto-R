##########################################################################################
#####################                                        			##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R             ##################
#####################                                       			##################
##########################################################################################


### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)

datos <- read.delim("C:/Users/ananc/Downloads/archivo-entrada.txt")

head(datos)
tail(datos)
dim(datos)

### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)

# Media y desviación estándar del Control: Las hacemos poniendo en los datos todas las réplicas correspondientes
#(columnas 2 y 4) aplicando la función fila a fila (1):

datos$Control_Mean <- apply(datos[,2:4],1,mean)
datos$Control_SD <- apply(datos[,2:4],1,sd)

# Media y desviación estándar de la Condición: Las hacemos poniendo en los datos todas las réplicas correspondientes
#(5 y 7):

datos$Condicion_Mean <- apply(datos[,5:7],1,mean)
datos$Condicion_SD <- apply(datos[,5:7],1,sd)

#Los resultados serán observados en el archivo: datos.

### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).

#Usamos par, que divide la gráfica en 1 fila y 2 columnas. 
#También es para que los dos histogramas salgan en una sola imagen:

par(mfrow=c(1,2))

#Para hacer el histograma Control:

hist(datos$Control_Mean, col="red", main="Media del control", xlab="Nivel de Expresión", ylab="Frecuencia" )
hist(datos$Condicion_Mean, col="blue", main="Media del Condición", xlab="Nivel de Expresión", ylab="Frecuencia" )

#Para volver a que se pongan las figuras de 1 en 1:

par(mfrow=c(1,1))

### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)

boxplot(datos$Control_Mean,datos$Condicion_Mean, names=c("Control","Condición"),
        col = c("red","blue"), main="Comparación boxplot: Control vs Condición", ylab="Expresión Media")

legend("topleft", legend=c("Control","Condicion"), fill=c("red","blue"))

### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla 

# Evaluación de la normalidad global con Shapiro-Wilk
shapiro.test(datos$Control_Mean)
shapiro.test(datos$Condicion_Mean)

# Aunque la distribución global no sigue una normalidad perfecta (p < 0.05),
# empleamos t.test de Student para comparar las 3 réplicas continuas de cada gen.
# No se usa wilcox.test porque con n=3 por grupo el p-valor mínimo alcanzable 
# matemáticamente es 0.1, imposibilitando detectar genes con p < 0.05.

p_valores <- apply(datos, 1, function(fila) {
  control_val   <- as.numeric(fila[2:4])
  condicion_val <- as.numeric(fila[5:7])
  
  # Si no hay variabilidad (sd = 0 en ambos), asignamos p-value = 1 para evitar error
  if (sd(control_val) == 0 && sd(condicion_val) == 0) {
    return(1)
  } else {
    res <- t.test(control_val, condicion_val)
    return(res$p.value)
  }
})

# Añadimos la columna de p-valores a la tabla
datos$p_value <- p_valores
head(datos)

### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)

#Primero, filtramos los genes significativos:

genes_diferenciales <- datos[datos$p_value<0.05,]

cat("Número total de genes diferencialmente expresados (pvalue <0.05) es:", nrow(genes_diferenciales))

#Ahora que sabemos que sí hay genes diferencialmente expresados, vemos cuáles se encuentran sobreexpresados
#y cuáles se encuentran reprimidos:

sobreexpresados <- genes_diferenciales[genes_diferenciales$Condicion_Mean > genes_diferenciales$Control_Mean, ]
reprimidos <- genes_diferenciales[genes_diferenciales$Condicion_Mean < genes_diferenciales$Control_Mean, ]

cat("Genes sobreexpresados:", nrow(sobreexpresados))
cat("Genes reprimidos:", nrow(reprimidos))

### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)

#Es necesario una matriz numérica. Se coge las 6 columnas de replicas y se
#ponen los nombres de gen como nombres de fila:

matriz_diferenciales <- as.matrix(genes_diferenciales[,2:7])
rownames(matriz_diferenciales) <- genes_diferenciales$Gene.Name

#Calculamos los dendogramas:

rowv = as.dendrogram(hclust(as.dist(1-cor(t(matriz_diferenciales)))))
colv = as.dendrogram(hclust(as.dist(1-cor(matriz_diferenciales))))

#Instalamos el paquete: gplots y, de ahí, obtenemos el heatmap:

library(gplots)
library(RColorBrewer)

heatmap.2(matriz_diferenciales, Rowv=rowv, Colv=colv, cexCol=0.7,
          col = rev(redblue(256)), scale = "row")
