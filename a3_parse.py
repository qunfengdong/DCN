nodes = dict()
edges = dict()

with open('all.edges.csv') as fh:
	next(fh)
	for line in fh:
		line = line.rstrip("\n|\r")
		cols = line.split(",")
		src = cols[0]
		target = cols[1]
		print(src, target)		
		if src not in nodes:
			nodes[src] = dict()
			nodes[src]['desc'] = cols[9].replace("'", "")
		if 'count' not in nodes[src]:
			nodes[src]['count'] = 0
		nodes[src]['count'] += 1
		if target not in nodes:
			nodes[target] = dict()
			nodes[target]['desc'] = cols[10].replace("'", "")
		if 'count' not in nodes[target]:
			nodes[target]['count'] = 0
		nodes[target]['count'] += 1
		edgeid = src + "_" + target
		if edgeid not in edges:
			edges[edgeid] = dict()
			edges[edgeid]['source'] = src
			edges[edgeid]['target'] = target
			edges[edgeid]['coef'] = cols[2]
			edges[edgeid]['exp_coef'] = cols[3]
			edges[edgeid]['se_coef'] = cols[4]
			edges[edgeid]['z'] = cols[5]
			edges[edgeid]['Pr'] = cols[6]
			edges[edgeid]['N'] = cols[7]
			edges[edgeid]['Pr_adjusted'] = cols[8]

datalist = []
for id in nodes:
	datalist.append("{{data:{{id:'{}',name:'{}',count:{}}}}}".format(id, nodes[id]['desc'], nodes[id]['count']))

for id in edges:
	datalist.append(" {{data:{{id:'{}',source:'{}',target:'{}',valueN:{},valuePrAdjusted:{},valueexp_coef:{}}}}}".format(id, edges[id]['source'], edges[id]['target'], edges[id]['N'],edges[id]['Pr_adjusted'],edges[id]['exp_coef']))

with open('all.edges.csv.js', 'w') as out:
	out.write("var jsonobj = [") 
	out.write(",".join(datalist))
	out.write("]")
