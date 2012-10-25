import re, fileinput, sys

def parse_line(content):
        try:
                ppe_re = re.search('H([A-Z])SDAFEQTSETIGV([A-Z])ANNA([A-Z])ND([A-Z])VRQRLL', content)
		ret_str = ""
		for i in range(1,5):
			ret_str = ret_str + ppe_re.group(i)
                return ret_str
        except:
                return None



for line in fileinput.input([sys.argv[1]]):
	ret = parse_line(line)
	if ret:
		print ret
