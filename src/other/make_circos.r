library(circlize)
x = rnorm(2300)
sectors = sample(letters[1:23], 2300, replace = TRUE)
circos.initialize(sectors, x=x)
circos.trackHist(sectors, x=x, col="#FF0000", bg.border="#FF0000", border="#FF0000", draw.density=TRUE)
circos.trackHist(sectors, x=x, col="#00FF00", bg.border="#00FF00", border="#00FF00")
circos.trackHist(sectors, x=x, col="#0000FF", bg.border="#0000FF", border="#0000FF", bin.size=0.5)
circos.trackHist(sectors, x=x, col="#C0C0C0", bg.border="#C0C0C0", border="#C0C0C0", bin.size=0.1)
