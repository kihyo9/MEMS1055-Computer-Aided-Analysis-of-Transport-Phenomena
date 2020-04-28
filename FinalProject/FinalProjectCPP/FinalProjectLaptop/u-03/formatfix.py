filename = input('Input file name (w/o ext): ')

f = open(filename + '.txt', 'r')
content = f.readlines()
f.close()

U = []
for line in content:
    data = list(map(float,line.split()))
    while(len(data) < 11):
        data.append(0.)
    U.append(data)

g = open(filename + "-fix.txt", 'w')
for line in U:
	first = True
	for el in line:
		if first:
			g.write(str(el))
			first = False
		else:
			g.write(" ")
			g.write(str(el))
	g.write("\n")
g.close()