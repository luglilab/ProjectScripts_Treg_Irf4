# libraries
import numpy as np
import matplotlib.pyplot as plt

# width of the bars
barWidth = 0.3

# Choose the height of the blue bars
bars1 = [5.29, 1.94, 52.54, 40.52]
# Choose the height of the cyan bars
bars2 = [18.66, 4.97, 31.36, 44.76]
#
bars3 = bars1 +bars2
#
label = ['5.29', '18.66', '1.94', '4.97', '52.54', '31.36', '40.52', '44.76']
# The x position of bars
r1 = np.arange(len(bars1))
r2 = [x + barWidth for x in r1]
r3 = r1- r2
# Create irf4 bars
plt.bar(r1, bars1, width=barWidth, color='orange', edgecolor='orange',  capsize=7, label='Irf4')
# Create batf bars
plt.bar(r2, bars2, width=barWidth, color='limegreen', edgecolor='limegreen',  capsize=7, label='Batf')

# general layout
plt.xticks([r + barWidth for r in range(len(bars1))], ['Promoter-TSS', 'Exon', 'Intron', 'Intergenic'])

plt.ylabel('Percentage(%)')
plt.legend()

plt.title("Peak Classification")
# Show graphic
plt.savefig("/Users/simonepuccio/Documents/HumanitasProjects/SP010/Figure4/batplot_annotation.eps", format='eps')
#plt.show()


import re
import periodic

symbol =[]
for i in range(len(periodic.elements)):
    symbol.append(periodic.elements[i].symbol.upper())

compounds = []
with open("teststefano.txt","r") as f:
    for line in f:
        compounds.append(line.strip("\n"))


for i in range(len(compounds)):
    for y in range(len(symbol)):
        if symbol[y] in compounds[i] and :


compounds_split = []
for i in range(len(compounds)):
    match = re.match(r"([a-z]+)([0-9]+)([a-z]+)([0-9]+)([a-z]+)([0-9]+)([a-z]+)([0-9]+)", compounds[i], re.I)
    if match:
        items = match.groups()
        compounds_split.append(items)



for i in range(len(compounds_split[0])):
    if "C" in compounds_split[0][i]:
        if (compounds_split[0][i + 1].isdigit() and compounds_split[0][i + 2].isdigit()):
            print compounds_split[0][i + 1] + compounds_split[0][i + 2]
        elif compounds_split[0][i + 1].isdigit():
            print compounds_split[0][i + 1]

string ="CH7O3"
for i in range(len(string)):
    if "C" in string[i]:
        if (string[i + 1].isdigit() and string[i + 2].isdigit()):
            print string[i + 1] + string[i + 2]
        elif string[i + 1].isdigit():
            print string[i + 1]
        else:
            print "1"