string = "matk, rbcL"
genes= string.split(",")

genes =[gene.replace(" ", "") for gene in genes]
print(genes)