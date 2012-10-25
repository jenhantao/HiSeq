import fileinput, sys

class Counter:
    def __init__(self):
        self.dict = {}
    def add(self, item):
        count = self.dict.get(item, 0)
        self.dict[item] = count + 1
    def counts(self, desc=None):
        """ Returns list of keys sorted by values.
        Pass desc as 1 if you want a descending sort. """
        result = map(None, self.dict.values(), self.dict.keys(  ))
        result.sort(  )
        if desc: result.reverse(  )
        return result

c = Counter( )
for line in fileinput.input([sys.argv[1]]):
	c.add(line)
#rrunning total
total = 0
for key in c.dict.keys():
	total = total + c.dict[key]
	#print key.rstrip()+"\t"+str(c.dict[key])
for key in c.dict.keys():
	print [sys.argv[1]],"\t",key.rstrip(),"\t", c.dict[key]/float(total)

