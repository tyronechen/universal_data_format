library(circlize)
x = rnorm(2300)
sectors = sample(letters[1:23], 2300, replace = TRUE)
circos.initialize(sectors, x=x)
circos.trackHist(sectors, x=x, col="#FF0000", border="#999999", draw.density=TRUE)
circos.trackHist(sectors, x=x, col="#00FF00", border="#999999")
circos.trackHist(sectors, x=x, col="#0000FF", border="#999999", bin.size=0.5)
circos.trackHist(sectors, x=x, col="#FF00FF", border="#999999", bin.size=0.1)
